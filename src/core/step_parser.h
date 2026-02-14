#pragma once

#include "types.h"
#include <string>
#include <TopoDS_Shape.hxx>

namespace camchamp {

class StepParser {
public:
    explicit StepParser(const std::string& step_file_path);
    ~StepParser();
    
    // Validate that STEP file exists and is valid format
    void Validate();
    
    // Detect units and confirm with user
    void DetectUnits();
    
    // Parse STEP file and load geometry
    void Parse();
    
    // Detect axis of revolution
    AxisOfRevolution DetectAxisOfRevolution();
    
    // Extract 2D profile from 3D geometry
    std::vector<Vertex> Extract2DProfile(const AxisOfRevolution& axis);
    
    // Get bounding box
    BoundingBox GetBoundingBox() const;
    
private:
    std::string step_file_path_;
    TopoDS_Shape shape_;
    BoundingBox bbox_;
};

bool IsParallel(const gp_Dir& dir1, const gp_Dir& dir2, double angular_tolerance_deg);
bool AxisPassesThroughPoint(const gp_Ax1& axis, const gp_Pnt& point, double distance_tolerance_mm);
bool AreSameAxis(const gp_Ax1& axis1, const gp_Ax1& axis2, double angular_tolerance_deg, double distance_tolerance_mm);

} // namespace camchamp
