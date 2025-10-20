#include "point_cloud.hpp"

void
visualize_point_cloud(const std::vector<Eigen::MatrixXd>& points_frames, const char* filename)
{
    // 点云数据
    // vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    // for ( const auto& frame : points_frames ) {
    //     for ( int i = 0; i < frame.rows(); ++i ) {
    //         points->InsertNextPoint(frame(i, 0), frame(i, 1), frame(i, 2));
    //     }
    // }
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for ( unsigned int i = 0; i < 10; ++i ) {
        points->InsertNextPoint(i, i, i);
    }
    // polyData对象
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    // 写入vtp文件
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(polyData);

    writer->Write();
}