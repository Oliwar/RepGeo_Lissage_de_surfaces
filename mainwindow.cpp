#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QVector3D>

// fonction pratique pour faire des tirages aléatoires
int randInt(int low, int high){return qrand() % ((high + 1) - low) + low;}

float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle(faceID);

    HalfedgeHandle heh = _mesh->halfedge_handle(fh);

    VertexHandle vh_A = _mesh->to_vertex_handle(heh);

    heh = _mesh->next_halfedge_handle(heh);
    VertexHandle vh_B = _mesh->to_vertex_handle(heh);

    heh = _mesh->next_halfedge_handle(heh);
    VertexHandle vh_C = _mesh->to_vertex_handle(heh);

    return ( (_mesh->point(vh_B) - _mesh->point(vh_A)) % (_mesh->point(vh_C) - _mesh->point(vh_A)) ).norm() / 2;
}

float MainWindow::areaBary(MyMesh* _mesh, int vertexID){
    VertexHandle vh = _mesh->vertex_handle(vertexID);

    float aireBary = 0;
    FaceHandle fh;
    for (MyMesh::VertexFaceIter vf_it=_mesh->vf_iter(vh); vf_it.is_valid(); ++vf_it)
    {
        fh = *vf_it;
        aireBary += faceArea(_mesh, fh.idx());
    }

    return aireBary /= 3.0f;
}

float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    FaceHandle fh = _mesh->face_handle(faceID);

    HalfedgeHandle heh = _mesh->halfedge_handle(fh);

    VertexHandle vhA;
    VertexHandle vhB;
    VertexHandle vhC;

    for(int i = 0; i < 3; ++i){
        vhA = _mesh->to_vertex_handle(heh);
        if(vhA == vh){
            vhB = _mesh->from_vertex_handle(heh);
            heh = _mesh->next_halfedge_handle(heh);
            vhC = _mesh->to_vertex_handle(heh);
            break;
        }
        heh = _mesh->next_halfedge_handle(heh);
    }

    return std::acos(
                (_mesh->point(vhB) - _mesh->point(vhA)).normalize()
                |
                (_mesh->point(vhC) - _mesh->point(vhA)).normalize()
                );
}

float MainWindow::cotan(float angle){
    return 1/(std::tan(angle));
}

Vec3f MainWindow::cotan_weight(MyMesh* _mesh, HalfedgeHandle heh){
    VertexHandle vh_i = _mesh->from_vertex_handle(heh);
    VertexHandle vh_j = _mesh->to_vertex_handle(heh);

    FaceHandle fh_b = _mesh->face_handle(heh);

    HalfedgeHandle heh_b = _mesh->prev_halfedge_handle(heh);
    VertexHandle vh_b = _mesh->from_vertex_handle(heh_b);

    heh = _mesh->opposite_halfedge_handle(heh);
    FaceHandle fh_a = _mesh->face_handle(heh);

    HalfedgeHandle heh_a = _mesh->next_halfedge_handle(heh);
    VertexHandle vh_a = _mesh->to_vertex_handle(heh_a);

    float cot_alpha = cotan(angleEE(_mesh, vh_a.idx(), fh_a.idx()));
    float cot_beta = cotan(angleEE(_mesh, vh_b.idx(), fh_b.idx()));
    float cot_alpha_beta = cot_alpha + cot_beta;

    Vec3f coord = _mesh->point(vh_j) - _mesh->point(vh_i);
    return cot_alpha_beta * coord;
}

Vec3f MainWindow::cotan_weight_sum(MyMesh* _mesh, int vertexID){
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    Vec3f total = Vec3f(0, 0, 0);
    for (MyMesh::VertexOHalfedgeIter voh_it=mesh.voh_iter(vh); voh_it.is_valid(); ++voh_it)
        total += cotan_weight(_mesh, *voh_it);
    return total;
}

std::vector<Vec3f> MainWindow::approximation_cotangentielle(MyMesh* _mesh){
    std::vector<Vec3f> tab(_mesh->n_vertices());
    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
        float area = 1/(2*areaBary(_mesh, v_it->idx()));
        Vec3f weight = cotan_weight_sum(_mesh, v_it->idx());
        Vec3f result = area * weight;
        tab.at(v_it->idx()) = result;
    }
    return tab;
}

bool MainWindow::are_neighbor(VertexHandle vh_i, VertexHandle vh_j){
    for (MyMesh::VertexVertexIter vv_it=mesh.vv_iter(vh_i); vv_it.is_valid(); ++vv_it)
        if(*vv_it == vh_j) return true;
    return false;
}

void MainWindow::matrice_Laplace_cotangentielle(MyMesh* _mesh){
    std::vector<std::vector<Vec3f>> mat_Laplace(_mesh->n_vertices(), std::vector<Vec3f>(_mesh->n_vertices()));
    VertexHandle vh_i;
    VertexHandle vh_j;
    for (MyMesh::VertexIter v_it=mesh.vertices_sbegin(); v_it!=mesh.vertices_end(); ++v_it) {
        vh_i = *v_it;
        for (MyMesh::VertexIter v_it2=mesh.vertices_sbegin(); v_it2!=mesh.vertices_end(); ++v_it2) {
            vh_j = *v_it2;
            if(vh_i == vh_j){
                mat_Laplace.at(vh_i.idx()).at(vh_j.idx()) = -cotan_weight_sum(_mesh, vh_i.idx());
            } else if(are_neighbor(vh_i, vh_j)){
                for(MyMesh::VertexOHalfedgeIter voh_it=mesh.voh_iter(vh_i); voh_it.is_valid(); ++voh_it){
                    if(mesh.to_vertex_handle(*voh_it) == vh_j){
                        mat_Laplace.at(vh_i.idx()).at(vh_j.idx()) = cotan_weight(_mesh, *voh_it);
                        break;
                    }
                }
            }
            else mat_Laplace.at(vh_i.idx()).at(vh_j.idx()) = Vec3f(0, 0, 0);
        }
    }
}

/* **** début de la partie boutons et IHM **** */

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

void MainWindow::on_pushButton_lissage_clicked()
{
    std::vector<Vec3f> tab = approximation_cotangentielle(&mesh);
    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
       mesh.point(mesh.vertex_handle(v_it->idx())) += h * lambda * tab[v_it->idx()];
    }
    displayMesh(&mesh);
}

Vec3f MainWindow::weight_sum(MyMesh* _mesh, int vertexID){
    VertexHandle vh = _mesh->vertex_handle(vertexID);
    Vec3f total = Vec3f(0, 0, 0);
    for (MyMesh::VertexVertexIter vv_it=mesh.vv_iter(vh); vv_it.is_valid(); ++vv_it){
        total += _mesh->point(*vv_it) - _mesh->point(vh);
    }
    return total;
}

std::vector<Vec3f> MainWindow::approximation_uniforme(MyMesh* _mesh){
    std::vector<Vec3f> tab(_mesh->n_vertices());
    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
        float area = 1/(2*areaBary(_mesh, v_it->idx()));
        Vec3f weight = weight_sum(_mesh, v_it->idx());
        Vec3f result = area * weight;
        tab.at(v_it->idx()) = result;
    }
    return tab;
}

void MainWindow::on_pushButton_lissage_uniforme_clicked()
{
    std::vector<Vec3f> tab = approximation_uniforme(&mesh);
    for (MyMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it){
       mesh.point(mesh.vertex_handle(v_it->idx())) += h * lambda * tab[v_it->idx()];
    }
    displayMesh(&mesh);
}
/* **** fin de la partie boutons et IHM **** */


/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_doubleSpinBox_h_valueChanged(double arg1)
{
    h = static_cast<float>(arg1);
}

void MainWindow::on_doubleSpinBox_lambda_valueChanged(double arg1)
{
    lambda = static_cast<float>(arg1);
}
