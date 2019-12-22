#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);

private slots:

    void on_pushButton_chargement_clicked();
    void on_pushButton_lissage_clicked();

    void on_doubleSpinBox_h_valueChanged(double arg1);

    void on_doubleSpinBox_lambda_valueChanged(double arg1);

    void on_pushButton_lissage_uniforme_clicked();

private:

    float faceArea(MyMesh* _mesh, int faceID);
    float areaBary(MyMesh* _mesh, int vertexID);
    float cotan(float angle);
    Vec3f cotan_weight(MyMesh *_mesh, HalfedgeHandle heh);
    Vec3f cotan_weight_sum(MyMesh* _mesh, int vertexID);
    std::vector<Vec3f> approximation_cotangentielle(MyMesh* _mesh);
    float angleEE(MyMesh* _mesh, int vertexID,  int faceID);
    bool are_neighbor(VertexHandle vh_i, VertexHandle vh_j);
    void matrice_Laplace_cotangentielle(MyMesh* _mesh);

    MyMesh mesh;

    float h = 0.01f;
    float lambda = 0.01f;

    Ui::MainWindow *ui;
    std::vector<Vec3f> approximation_uniforme(MyMesh *_mesh);
    Vec3f weight_sum(MyMesh *_mesh, int vertexID);
};

#endif // MAINWINDOW_H
