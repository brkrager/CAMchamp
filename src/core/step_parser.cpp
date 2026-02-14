#include "step_parser.h"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>

//needed to open and validate step file
#include <STEPControl_Controller.hxx>
#include <STEPControl_Reader.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>

//needed to get bounding box
#include <TopoDS_Shape.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>

//needed to identify surface types
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

//needed to detect revolution axis
#include <gp_Cylinder.hxx>
#include <gp_Cone.hxx>
#include <gp_Torus.hxx>
#include <gp_Ax1.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Sphere.hxx>

namespace camchamp {

    struct CandidateAxis {
        gp_Ax1 axis;
        double weight;
        std::vector<std::string> contributing_surfaces; // For debugging

        CandidateAxis(const gp_Ax1& ax, double w, const std::string& type)
            : axis(ax), weight(w) {
            contributing_surfaces.push_back(type);
        }

        void addWeight(double w, const std::string& type) {
            weight += w;
            contributing_surfaces.push_back(type);
        }
    };
    struct PlanarSurface {
        gp_Dir normal;
        double area;

        PlanarSurface(const gp_Dir& n, double a) : normal(n), area(a) {}
    };
    struct SphericalSurface {
        gp_Pnt center;
        double area;

        SphericalSurface(const gp_Pnt& c, double a) : center(c), area(a) {}
    };
    struct GroupedAxis {
        gp_Ax1 representative_axis;  // Use first axis in group as representative
        double total_weight;
        std::vector<std::string> all_contributing_surfaces;
        int num_axes_in_group;

        GroupedAxis(const gp_Ax1& ax, double w, const std::vector<std::string>& surfaces)
            : representative_axis(ax), total_weight(w),
            all_contributing_surfaces(surfaces), num_axes_in_group(1) {
        }

        void mergeAxis(const CandidateAxis& other) {
            total_weight += other.weight;
            all_contributing_surfaces.insert(all_contributing_surfaces.end(),
                other.contributing_surfaces.begin(),
                other.contributing_surfaces.end());
            num_axes_in_group++;
        }
    };

StepParser::StepParser(const std::string& step_file_path)
    : step_file_path_(step_file_path)
{
}

StepParser::~StepParser() {
}

void StepParser::Validate() {
    std::cout << "  Validating STEP file: " << step_file_path_ << "\n";

    std::ifstream test(step_file_path_);
    if (!test.good()) {
        throw std::runtime_error("  File doesn't exist or can't be opened!");
    }
    test.close();
}

void StepParser::DetectUnits() {
    std::cout << "  Detecting units from STEP file...\n";

    // Verify we have a valid shape to work with
    if (shape_.IsNull()) {
        throw std::runtime_error("  No geometry loaded. Call Parse() before DetectUnits().");
    }

    // Calculate bounding box using OpenCASCADE
    Bnd_Box bnd_box;
    BRepBndLib::Add(shape_, bnd_box);

    // Check if bounding box is valid
    if (bnd_box.IsVoid()) {
        throw std::runtime_error("  Failed to calculate bounding box - shape may be empty.");
    }

    // Extract min/max coordinates
    double xmin, ymin, zmin, xmax, ymax, zmax;
    bnd_box.Get(xmin, ymin, zmin, xmax, ymax, zmax);

    // Store in member variable
    bbox_.min_x = xmin;
    bbox_.max_x = xmax;
    bbox_.min_y = ymin;
    bbox_.max_y = ymax;
    bbox_.min_z = zmin;
    bbox_.max_z = zmax;
  
    std::cout << "  Detected dimensions: "
        << bbox_.GetWidth() << "mm x "
        << bbox_.GetDepth() << "mm x "
        << bbox_.GetHeight() << "mm\n";
}

void StepParser::Parse() {
    std::cout << "  Parsing STEP geometry...\n";
    
    STEPControl_Controller::Init();
    STEPControl_Reader reader;
    IFSelect_ReturnStatus status = reader.ReadFile(step_file_path_.c_str());

    if (status != IFSelect_RetDone) {
        throw std::runtime_error("  Failed to read STEP file: Invalid format, exit code " + std::to_string(status));
    }
    reader.TransferRoots();
    shape_ = reader.OneShape();
    std::cout << "  File parsed\n";

    // Count solid bodies
    int solid_count = 0;
    TopExp_Explorer explorer(shape_, TopAbs_SOLID);
    for (; explorer.More(); explorer.Next()) {
        solid_count++;
    }
    if (solid_count < 1) {
        throw std::runtime_error("  No solid bodies found in file");
    }
    else if (solid_count > 1) {
        throw std::runtime_error("  Found " + std::to_string(solid_count) + " solid bodies. Upload only 1 solid body per file.");
    }
    

    /*TColStd_SequenceOfAsciiString lengthUnits;
    TColStd_SequenceOfAsciiString angleUnits;
    TColStd_SequenceOfAsciiString solidAngleUnits;
    reader.FileUnits(lengthUnits, angleUnits, solidAngleUnits);

    for (int i = 1; i <= lengthUnits.Length(); ++i) {
        std::cout << "  Length unit found: " << lengthUnits.Value(i).ToCString() << "\n";
    }
    for (int i = 1; i <= angleUnits.Length(); ++i) {
        std::cout << "  Angle unit found: " << angleUnits.Value(i).ToCString() << "\n";
    }
    for (int i = 1; i <= solidAngleUnits.Length(); ++i) {
        std::cout << "  Solid angle unit found: " << solidAngleUnits.Value(i).ToCString() << "\n";
    }*/

}

AxisOfRevolution StepParser::DetectAxisOfRevolution() {

    const double ANGULAR_TOLERANCE = 0.5;  // degrees
    const double DISTANCE_TOLERANCE = 0.02; // mm
    const double WEIGHT_THRESHOLD_PERCENT = 0.0; //the winning axis must have at least this much % of total surface area to continue

    // Storage for surfaces
    std::vector<CandidateAxis> candidate_axes;
    std::vector<PlanarSurface> planar_surfaces;
    std::vector<SphericalSurface> spherical_surfaces;
    
    // Counters for each surface type
    int cylindrical_count = 0;
    int conical_count = 0;
    int planar_count = 0;
    int spherical_count = 0;
    int toroidal_count = 0;
    int other_count = 0;

    gp_Ax1 axis;
    gp_Pnt loc;
    gp_Dir dir;
    gp_Cylinder cylinder;
    gp_Cone cone;
    gp_Torus torus;
    gp_Pln plane;
    gp_Sphere sphere;
    double totalSurfaceArea = 0.0;

    // Iterate through all faces in the shape
    TopExp_Explorer face_explorer(shape_, TopAbs_FACE);
    for (; face_explorer.More(); face_explorer.Next()) {
        TopoDS_Face face = TopoDS::Face(face_explorer.Current());

        // Adapt the face to query its geometry
        BRepAdaptor_Surface surface(face);
        GeomAbs_SurfaceType surf_type = surface.GetType();

        // Calculate surface area
        GProp_GProps props;
        BRepGProp::SurfaceProperties(face, props);
        double area = props.Mass();
        totalSurfaceArea += area;

        // Classify and count
        switch (surf_type) {
        case GeomAbs_Cylinder:
            cylindrical_count++;
            cylinder = surface.Cylinder();
            axis = cylinder.Axis();

            loc = axis.Location();
            dir = axis.Direction();

            //std::cout << "    Cylindrical surface (area: " << area << ")\n";
            //std::cout << "      Axis location: (" << loc.X() << ", " << loc.Y() << ", " << loc.Z() << ")\n";
            //std::cout << "      Axis direction: (" << dir.X() << ", " << dir.Y() << ", " << dir.Z() << ")\n";

            candidate_axes.emplace_back(axis, area, "Cylinder");
            break;
        case GeomAbs_Cone:
            conical_count++;
            cone = surface.Cone();
            axis = cone.Axis();

            loc = axis.Location();
            dir = axis.Direction();

            //std::cout << "    Conical surface (area: " << area << ")\n";
            //std::cout << "      Axis location: (" << loc.X() << ", " << loc.Y() << ", " << loc.Z() << ")\n";
            //std::cout << "      Axis direction: (" << dir.X() << ", " << dir.Y() << ", " << dir.Z() << ")\n";

            candidate_axes.emplace_back(axis, area, "Cone");
            break;
        case GeomAbs_Torus:
            toroidal_count++;
            torus = surface.Torus();
            axis = torus.Axis();

            loc = axis.Location();
            dir = axis.Direction();

            //std::cout << "    Toroidal surface (area: " << area << ")\n";
            //std::cout << "      Axis location: (" << loc.X() << ", " << loc.Y() << ", " << loc.Z() << ")\n";
            //std::cout << "      Axis direction: (" << dir.X() << ", " << dir.Y() << ", " << dir.Z() << ")\n";

            candidate_axes.emplace_back(axis, area, "Torus");
            break;
        case GeomAbs_Plane:
            planar_count++;
            plane = surface.Plane();
            dir = plane.Axis().Direction();

            //std::cout << "    Planar surface (area: " << area << ")\n";
            //std::cout << "      Normal: (" << dir.X() << ", " << dir.Y() << ", " << dir.Z() << ")\n";

            planar_surfaces.emplace_back(dir, area);
            break;
        case GeomAbs_Sphere:
            spherical_count++;
            sphere = surface.Sphere();
            loc = sphere.Location();

            //std::cout << "    Spherical surface (area: " << area << ")\n";
            //std::cout << "      Center: (" << loc.X() << ", " << loc.Y() << ", " << loc.Z() << ")\n";

            spherical_surfaces.emplace_back(loc, area);
            break;
        default:
            other_count++;
            //std::cout << "    Other surface type: " << surf_type << " (area: " << area << ")\n";
            break;
        }
    }

    // Add weights from planar surfaces
    for (const auto& plane : planar_surfaces) {
        //std::cout << "    Checking planar surface (area: " << plane.area << ")...\n";
        bool matched = false;

        for (auto& candidate : candidate_axes) {
            if (IsParallel(plane.normal, candidate.axis.Direction(), ANGULAR_TOLERANCE)) {
                //std::cout << "      Matched to candidate axis (parallel normal)\n";
                candidate.addWeight(plane.area, "Plane");
                matched = true;
                break; // Only add to first matching axis
            }
        }

        if (!matched) {
            //std::cout << "      No matching axis found\n";
        }
    }

    // Add weights from spherical surfaces
    for (const auto& sphere : spherical_surfaces) {
        //std::cout << "    Checking spherical surface (area: " << sphere.area << ")...\n";
        bool matched = false;

        for (auto& candidate : candidate_axes) {
            if (AxisPassesThroughPoint(candidate.axis, sphere.center, DISTANCE_TOLERANCE)) {
                //std::cout << "      Matched to candidate axis (axis passes through center)\n";
                candidate.addWeight(sphere.area, "Sphere");
                matched = true;
                break; // Only add to first matching axis
            }
        }

        if (!matched) {
            //std::cout << "      No matching axis found\n";
        }
    }


        // Summary
    /*std::cout << "  Surface summary:\n";
    std::cout << "    Cylindrical: " << cylindrical_count << "\n";
    std::cout << "    Conical: " << conical_count << "\n";
    std::cout << "    Planar: " << planar_count << "\n";
    std::cout << "    Spherical: " << spherical_count << "\n";
    std::cout << "    Toroidal: " << toroidal_count << "\n";
    std::cout << "    Other: " << other_count << "\n";*/



    //  Group candidate axes that are the same
    std::vector<GroupedAxis> grouped_axes;

    for (const auto& candidate : candidate_axes) {
        bool found_group = false;

        // Check if this candidate matches any existing group
        for (auto& group : grouped_axes) {
            if (AreSameAxis(candidate.axis, group.representative_axis,
                ANGULAR_TOLERANCE, DISTANCE_TOLERANCE)) {
                //std::cout << "    Merging axis into existing group (weight: "
                //    << candidate.weight << ")\n";
                group.mergeAxis(candidate);
                found_group = true;
                break;
            }
        }

        // If no matching group found, create a new group
        if (!found_group) {
            //std::cout << "    Creating new group (weight: " << candidate.weight << ")\n";
            grouped_axes.emplace_back(candidate.axis, candidate.weight, candidate.contributing_surfaces);
        }
    }

    //display all unique axis options
    //std::cout << "\n  Grouped axes:\n";
    for (size_t i = 0; i < grouped_axes.size(); i++) {
        const auto& group = grouped_axes[i];
        loc = group.representative_axis.Location();
        dir = group.representative_axis.Direction();

        /*std::cout << "    Group " << i + 1 << ":\n";
        std::cout << "      Number of axes merged: " << group.num_axes_in_group << "\n";
        std::cout << "      Representative axis location: (" << loc.X() << ", "
            << loc.Y() << ", " << loc.Z() << ")\n";
        std::cout << "      Representative axis direction: (" << dir.X() << ", "
            << dir.Y() << ", " << dir.Z() << ")\n";
        std::cout << "      Total combined weight: " << group.total_weight << "\n";*/
    }

    // Find the axis with maximum weight
    if (grouped_axes.empty()) {
        throw std::runtime_error("No candidate axes found (no cylindrical, conical, or toroidal surfaces)");
    }

    size_t winning_index = 0;
    double max_weight = grouped_axes[0].total_weight;

    for (size_t i = 1; i < grouped_axes.size(); i++) {
        if (grouped_axes[i].total_weight > max_weight) {
            max_weight = grouped_axes[i].total_weight;
            winning_index = i;
        }
    }

    const GroupedAxis& winner = grouped_axes[winning_index];
    double weight_percentage = (max_weight / totalSurfaceArea) * 100.0;

    //std::cout << "    Winning axis: Group " << (winning_index + 1) << "\n";
    //std::cout << "    Weight: " << max_weight << " (" << weight_percentage << "% of total)\n";

    if (weight_percentage < WEIGHT_THRESHOLD_PERCENT) {
        throw std::runtime_error("No clear axis of revolution detected. "
            "Winning axis only accounts for " +
            std::to_string(weight_percentage) +
            "% of surface area (threshold: " +
            std::to_string(WEIGHT_THRESHOLD_PERCENT) + "%)");
    }

    //std::cout << "    Validation passed (exceeds " << WEIGHT_THRESHOLD_PERCENT << "% threshold)\n";
    
    // Create and return the winning axis
    AxisOfRevolution result;
    dir = winner.representative_axis.Direction();
    loc = winner.representative_axis.Location();

    result.direction_x = dir.X();
    result.direction_y = dir.Y();
    result.direction_z = dir.Z();
    result.location_x = loc.X();
    result.location_y = loc.Y();
    result.location_z = loc.Z();
    result.total_weight = winner.total_weight;

    std::cout << "  Direction: (" << result.direction_x << ", "
        << result.direction_y << ", " << result.direction_z << ")\n";
    std::cout << "  Location: (" << result.location_x << ", "
        << result.location_y << ", " << result.location_z << ")\n";
    std::cout << "  Represented surface area: " << weight_percentage << "%\n";

    return result;
}



std::vector<Vertex> StepParser::Extract2DProfile(const AxisOfRevolution& axis) {
    std::cout << "Extracting 2D profile...\n";
    
    // TODO: Implement profile extraction
    
    // Placeholder: Simple shaft profile
    std::vector<Vertex> profile;
    profile.push_back(Vertex(4.0, 0.5));    // Front face, OD
    profile.push_back(Vertex(1.25, 0.5));   // End of OD
    profile.push_back(Vertex(1.25, 0.25));  // Start of ID
    profile.push_back(Vertex(4.0, 0.25));   // End of ID (back to front)
    
    std::cout << "  [PLACEHOLDER] Extracted " << profile.size() << " vertices\n";
    return profile;
}

BoundingBox StepParser::GetBoundingBox() const {
    return bbox_;
}

// Check if a direction is parallel to an axis direction (within angular tolerance)
bool IsParallel(const gp_Dir& dir1, const gp_Dir& dir2, double angular_tolerance_deg = 0.5) {
    double angle_rad = dir1.Angle(dir2);
    double angle_deg = angle_rad * 180.0 / M_PI;

    // Check if parallel (0 degrees) or anti-parallel (180 degrees)
    return (angle_deg < angular_tolerance_deg) || (angle_deg > (180.0 - angular_tolerance_deg));
}

// Check if an axis passes through a point (within distance tolerance)
bool AxisPassesThroughPoint(const gp_Ax1& axis, const gp_Pnt& point, double distance_tolerance_mm = 0.02) {
    gp_Lin line(axis);
    double distance = line.Distance(point);
    return distance < distance_tolerance_mm;
}

// Check if two axes are the same (parallel directions and collinear)
bool AreSameAxis(const gp_Ax1& axis1, const gp_Ax1& axis2,
    double angular_tolerance_deg = 0.5,
    double distance_tolerance_mm = 0.02) {
    // Check if directions are parallel
    if (!IsParallel(axis1.Direction(), axis2.Direction(), angular_tolerance_deg)) {
        return false;
    }

    // Check if axes are collinear (axis1 passes through a point on axis2)
    gp_Pnt point_on_axis2 = axis2.Location();
    return AxisPassesThroughPoint(axis1, point_on_axis2, distance_tolerance_mm);
}

} // namespace camchamp



