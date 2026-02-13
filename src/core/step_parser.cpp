#include "step_parser.h"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>

#include <STEPControl_Controller.hxx>
#include <STEPControl_Reader.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

namespace camchamp {

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
    std::cout << "  Detecting axis of revolution...\n";
    
    // Counters for each surface type
    int cylindrical_count = 0;
    int conical_count = 0;
    int planar_count = 0;
    int spherical_count = 0;
    int toroidal_count = 0;
    int other_count = 0;

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

        // Classify and count
        switch (surf_type) {
        case GeomAbs_Cylinder:
            cylindrical_count++;
            std::cout << "    Cylindrical surface (area: " << area << ")\n";
            break;
        case GeomAbs_Cone:
            conical_count++;
            std::cout << "    Conical surface (area: " << area << ")\n";
            break;
        case GeomAbs_Plane:
            planar_count++;
            std::cout << "    Planar surface (area: " << area << ")\n";
            break;
        case GeomAbs_Sphere:
            spherical_count++;
            std::cout << "    Spherical surface (area: " << area << ")\n";
            break;
        case GeomAbs_Torus:
            toroidal_count++;
            std::cout << "    Toroidal surface (area: " << area << ")\n";
            break;
        default:
            other_count++;
            std::cout << "    Other surface type: " << surf_type << " (area: " << area << ")\n";
            break;
        }
    }

    // Summary
    std::cout << "  Surface summary:\n";
    std::cout << "    Cylindrical: " << cylindrical_count << "\n";
    std::cout << "    Conical: " << conical_count << "\n";
    std::cout << "    Planar: " << planar_count << "\n";
    std::cout << "    Spherical: " << spherical_count << "\n";
    std::cout << "    Toroidal: " << toroidal_count << "\n";
    std::cout << "    Other: " << other_count << "\n";


    
    AxisOfRevolution axis;
    axis.direction_x = 0;
    axis.direction_y = 0;
    axis.direction_z = 1;  // Z-axis
    axis.location_x = 0;
    axis.location_y = 0;
    axis.location_z = 0;
    axis.total_weight = 100.0;  // Placeholder
    
    std::cout << "  [PLACEHOLDER] Axis: Z-axis (0, 0, 1)\n";
    return axis;
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

} // namespace camchamp
