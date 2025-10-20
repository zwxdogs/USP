#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <eigen3/Eigen/Dense>

void
visualize_point_cloud(const std::vector<Eigen::MatrixXd>& points_frames, const char* filename);