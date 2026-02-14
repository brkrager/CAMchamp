#pragma once

#include <string>
#include <vector>
#include <memory>

namespace camchamp {

// Forward declarations
struct Vertex;
struct Feature;
struct Operation;
struct Tool;
struct Stock;
struct Material;
struct MachineConfig;

// ============================================================================
// ENUMERATIONS
// ============================================================================

enum class MaterialFamily {
    Unknown,
    Plastic,
    Aluminum,
    Brass,
    Steel_Mild,
    Steel_Hardened,
    Stainless_Austenitic,
    Stainless_Ferritic,
    Cast_Iron,
    Titanium,
    Inconel
};

enum class OperationType {
    Face,
    CenterDrill,
    Drill,
    RoughTurnOD,
    FinishTurnOD,
    RoughBoreID,
    FinishBoreID,
    GrooveOD,
    GrooveID,
    ThreadOD,
    ThreadID,
    PartOff,
    Unknown
};

enum class FeatureType {
    Face,
    OD_Straight,
    OD_Taper,
    OD_Radius,
    OD_Chamfer,
    ID_Straight,
    ID_Taper,
    ID_Radius,
    Groove,
    Thread,
    Unknown
};

enum class ToolType {
    TurningInsert,
    BoringBar,
    Drill,
    CenterDrill,
    PartingTool,
    GroovingTool,
    ThreadingTool,
    Endmill,
    Unknown
};

enum class ToolUsage {
    RoughOnly,
    FinishOnly,
    Both
};

enum class SeverityLevel {
    Error,
    Warning,
    Info,
    Debug
};

// ============================================================================
// BASIC STRUCTURES
// ============================================================================

struct Vertex {
    double x;  // Z-position in lathe coordinates (inches)
    double y;  // Radius (inches)
    
    Vertex() : x(0.0), y(0.0) {}
    Vertex(double x_, double y_) : x(x_), y(y_) {}
};

struct BoundingBox {
    double min_x, max_x;
    double min_y, max_y;
    double min_z, max_z;
    
    BoundingBox() 
        : min_x(0), max_x(0)
        , min_y(0), max_y(0)
        , min_z(0), max_z(0) {}
    
    double GetWidth() const { return max_x - min_x; }
    double GetDepth() const { return max_y - min_y; }
    double GetHeight() const { return max_z - min_z; }
};

struct AxisOfRevolution {
    double direction_x, direction_y, direction_z;  // Unit vector
    double location_x, location_y, location_z;     // Point on axis
    double total_weight;                           // Total surface area
    
    AxisOfRevolution()
        : direction_x(0), direction_y(0), direction_z(1)
        , location_x(0), location_y(0), location_z(0)
        , total_weight(0) {}
};



// ============================================================================
// MATERIAL
// ============================================================================

struct Material {
    std::string material_id;
    std::string description;
    MaterialFamily material_family;
    
    // Physical properties
    double density_lb_per_in3;
    double elastic_modulus_psi;
    double shear_modulus_psi;
    double poissons_ratio;
    double tensile_strength_psi;
    double yield_strength_psi;
    double shear_strength_psi;
    double hardness_brinell;
    
    // Machining properties
    double machinability_rating;
    double specific_cutting_force_psi;
    
    // Stock safety factor
    double stock_diameter_safety_factor;
    
    Material()
        : material_family(MaterialFamily::Unknown)
        , density_lb_per_in3(0)
        , elastic_modulus_psi(0)
        , shear_modulus_psi(0)
        , poissons_ratio(0)
        , tensile_strength_psi(0)
        , yield_strength_psi(0)
        , shear_strength_psi(0)
        , hardness_brinell(0)
        , machinability_rating(0)
        , specific_cutting_force_psi(0)
        , stock_diameter_safety_factor(1.03) {}
};

// ============================================================================
// STOCK
// ============================================================================

struct Stock {
    std::string stock_id;
    std::string material_id;
    double diameter_in;
    double length_in;
    double stickout_in;
    bool tailstock_engaged;
    
    Stock()
        : diameter_in(0)
        , length_in(0)
        , stickout_in(0)
        , tailstock_engaged(false) {}
};

// ============================================================================
// TOOL
// ============================================================================

struct ToolMaterialParams {
    std::string material_family;
    double cutting_speed_sfm;
    double feed_per_rev_in;
    double max_depth_of_cut_in;
    double chipload_min_in;
    double chipload_max_in;
    std::string insert_grade;
    
    ToolMaterialParams()
        : cutting_speed_sfm(0)
        , feed_per_rev_in(0)
        , max_depth_of_cut_in(0)
        , chipload_min_in(0)
        , chipload_max_in(0) {}
};

struct Tool {
    std::string tool_id;
    std::string description;
    ToolType tool_type;
    int tool_number;
    
    // Geometry
    double nose_radius_in;
    double cutting_edge_angle_deg;
    
    // Holder
    std::string holder_type;
    double shank_width_in;
    double shank_height_in;
    double overhang_length_in;
    double cutting_edge_height_in;
    
    // Rigidity
    double tool_stiffness_N_per_mm;
    double max_cutting_force_N;
    double deflection_limit_in;
    
    // Usage
    ToolUsage usage;
    std::vector<std::string> material_whitelist;
    std::vector<std::string> material_blacklist;
    
    // Material-specific parameters
    std::vector<ToolMaterialParams> material_params;
    
    Tool()
        : tool_type(ToolType::Unknown)
        , tool_number(0)
        , nose_radius_in(0)
        , cutting_edge_angle_deg(0)
        , shank_width_in(0)
        , shank_height_in(0)
        , overhang_length_in(0)
        , cutting_edge_height_in(0)
        , tool_stiffness_N_per_mm(0)
        , max_cutting_force_N(0)
        , deflection_limit_in(0.002)
        , usage(ToolUsage::Both) {}
    
    // Helper methods
    ToolMaterialParams GetParamsForMaterial(const std::string& material_family) const;
    bool IsCompatibleWithMaterial(const std::string& material_family) const;
};

// ============================================================================
// FEATURE
// ============================================================================

struct Feature {
    FeatureType type;
    std::vector<Vertex> profile;
    double Z_start_in;
    double Z_end_in;
    double diameter_start_in;
    double diameter_end_in;
    bool accessible_in_setup;
    
    Feature()
        : type(FeatureType::Unknown)
        , Z_start_in(0)
        , Z_end_in(0)
        , diameter_start_in(0)
        , diameter_end_in(0)
        , accessible_in_setup(true) {}
};

// ============================================================================
// OPERATION
// ============================================================================

struct CuttingParameters {
    double cutting_speed_sfm;
    double spindle_rpm;
    double feed_per_rev_in;
    double feed_rate_ipm;
    double depth_of_cut_in;
    int number_of_passes;
    double stock_to_leave_in;
    double finish_stock_to_leave_in;
    
    CuttingParameters()
        : cutting_speed_sfm(0)
        , spindle_rpm(0)
        , feed_per_rev_in(0)
        , feed_rate_ipm(0)
        , depth_of_cut_in(0)
        , number_of_passes(1)
        , stock_to_leave_in(0)
        , finish_stock_to_leave_in(0) {}
};

struct CalculatedValues {
    double material_removal_volume_in3;
    double estimated_cycle_time_sec;
    double max_cutting_force_lbf;
    double calculated_deflection_in;
    double horsepower_required;
    double metal_removal_rate_in3_per_min;
    double surface_finish_ra_microinch;
    
    CalculatedValues()
        : material_removal_volume_in3(0)
        , estimated_cycle_time_sec(0)
        , max_cutting_force_lbf(0)
        , calculated_deflection_in(0)
        , horsepower_required(0)
        , metal_removal_rate_in3_per_min(0)
        , surface_finish_ra_microinch(0) {}
};

struct Operation {
    std::string operation_id;
    OperationType operation_type;
    std::string operation_name;
    int sequence_number;
    int setup_number;
    
    std::shared_ptr<Tool> tool;
    std::string spindle;
    
    CuttingParameters parameters;
    CalculatedValues calculated_values;
    
    // Geometry (varies by operation type)
    double start_diameter_in;
    double end_diameter_in;
    double Z_start_in;
    double Z_end_in;
    
    std::string coolant;
    std::vector<std::string> warnings;
    std::string comment;
    
    Operation()
        : operation_type(OperationType::Unknown)
        , sequence_number(0)
        , setup_number(1)
        , start_diameter_in(0)
        , end_diameter_in(0)
        , Z_start_in(0)
        , Z_end_in(0)
        , coolant("flood") {}
};

// ============================================================================
// MACHINE CONFIGURATION
// ============================================================================

struct AxisConfig {
    double travel_in;
    double max_feed_ipm;
    double rapid_ipm;
    
    AxisConfig()
        : travel_in(0), max_feed_ipm(0), rapid_ipm(0) {}
};

struct SpindleConfig {
    std::string type;  // "lathe" or "mill"
    double max_rpm;
    double min_rpm;
    double max_power_hp;
    double max_torque_ftlb;
    
    SpindleConfig()
        : max_rpm(0), min_rpm(0), max_power_hp(0), max_torque_ftlb(0) {}
};

struct PartOffConfig {
    bool part_off_break;
    double part_off_break_radius_in;
    bool part_off_auto_annular_ring;
    double part_off_annular_ring_thickness_in;
    double spindle_deceleration_rad_per_sec2;
    double annular_ring_safety_factor;
    
    PartOffConfig()
        : part_off_break(false)
        , part_off_break_radius_in(-0.025)
        , part_off_auto_annular_ring(false)
        , part_off_annular_ring_thickness_in(0.030)
        , spindle_deceleration_rad_per_sec2(10.0)
        , annular_ring_safety_factor(3.0) {}
};

struct MachineConfig {
    std::string machine_id;
    std::string machine_type;
    std::string description;
    
    AxisConfig X_axis;
    AxisConfig Z_axis;
    SpindleConfig main_spindle;
    
    bool tailstock_available;
    double tailstock_max_force_lbf;
    
    bool part_catcher_available;
    double part_catcher_max_length_in;
    double part_catcher_max_diameter_in;
    
    PartOffConfig part_off_config;
    
    MachineConfig()
        : tailstock_available(false)
        , tailstock_max_force_lbf(0)
        , part_catcher_available(false)
        , part_catcher_max_length_in(0)
        , part_catcher_max_diameter_in(0) {}
};

// ============================================================================
// PROGRAM OUTPUT
// ============================================================================

struct ProgramMetadata {
    std::string part_name;
    std::string step_file_path;
    std::string material;
    std::string machine;
    Stock stock;
    std::string generated_by;
    std::string generated_date;
    int total_operations;
    double estimated_total_cycle_time_min;
};

struct ProgramOutput {
    ProgramMetadata metadata;
    std::vector<Operation> operations;
    std::vector<std::string> warnings;
    
    // Summary statistics
    double total_material_removed_in3;
    double total_cycle_time_min;
    int total_tool_changes;
    double max_horsepower_required;
    double stock_utilization_percent;
};

} // namespace camchamp
