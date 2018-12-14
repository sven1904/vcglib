/****************************************************************************
* (C) Sven Tittel, 2019                                                     *
*                                                                           *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/*! \file urdf_from_file.cpp
\ingroup code_sample

\brief An example of computing the inertia properties of meshes

Two meshes are created a rectangular box and a torus and their mass properties are computed and shown.
The result should match the closed formula for these objects (with a reasonable approximation)

*/

#include <fstream>
#include <vcg/complex/complex.h>
#include <GL/gl.h>
#include <wrap/io_trimesh/import_dae.h>
#include <wrap/io_trimesh/import.h>
#include <vcg/complex/algorithms/inertia.h>
#include <vcg/complex/algorithms/create/platonic.h>

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>   ::AsVertexType,
                                            vcg::Use<MyEdge>     ::AsEdgeType,
                                            vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags> {};
class MyFace    : public vcg::Face<MyUsedTypes, vcg::face::FFAdj, vcg::face::Normal3f, vcg::face::VertexRef, vcg::face::BitFlags> {};
class MyEdge    : public vcg::Edge<MyUsedTypes> {};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace>, std::vector<MyEdge>> {};

int main( int argc, char **argv )
{
  std::vector<double> VolPerLink;
  std::vector<vcg::Matrix33f> ITvector;
  std::vector<vcg::Point3f> COMvector;
  std::vector<Eigen::VectorXf> joints;
  vcg::tri::io::InfoDAE info;
  double mass = 1.0;
  double volume = 0.0;

  if (argc <= 1)
  {
      printf("No mesh file (.dae) provided!\n");
      return -1;
  }


  for (int i = 1; i < argc; i++)
  {
      double tmp = atof(argv[i]);
      if (tmp > 0)
      {
          mass = tmp;
          printf("Overall mass is: %f kg\n", mass);
          continue;
      }

      size_t len = strlen(argv[i]);
      if (len > 3 && strcmp(argv[i]+len-4, ".txt") == 0)
      {
          joints.clear();
          printf("Read file '%s' as joint transformation info.\n", argv[i]);
          std::ifstream file(argv[i]);
          char buffer[32];
          char* value;
          while (file.getline(buffer, 32, '\n'))
          {
              printf("char buffer: %s\n", buffer);
              Eigen::VectorXf joint;
              joint.setZero(6);
              value = std::strtok(buffer, " ");
              int index = 0;
              while (value)
              {
                  joint[index++] = static_cast<float>(atof(value));
                  value = std::strtok(nullptr, " ");
              }
              joints.push_back(joint);
          }
          file.close();
          continue;
      }

      MyMesh linkMesh;
      if (len > 3 && strcmp(argv[i]+len-4, ".dae") == 0)
      {
          if (vcg::tri::io::ImporterDAE<MyMesh>::Open(linkMesh, argv[i], info) != 0)
          {
              printf("Could not open file '%s'\n", argv[i]);
              return -1;
          }
      }
      else
      {
          if (vcg::tri::io::Importer<MyMesh>::Open(linkMesh, argv[i]) != 0)
          {
              printf("Could not open file '%s'\n", argv[i]);
              return -1;
          }
      }

      vcg::tri::Inertia<MyMesh> Ib(linkMesh);
      COMvector.push_back(Ib.CenterOfMass());
      vcg::Matrix33f ITensor;

      double vol = static_cast<double>(std::abs(Ib.Mass()));
      VolPerLink.push_back(vol);
      printf("Volume: %14.11f + %14.11f", volume, vol);
      volume += vol;
      printf(" = %14.11f\n", volume);
      Ib.InertiaTensor(ITensor);
      ITvector.push_back(ITensor);
  }

  printf("URDF data for %d links with overall mass of %.3f kg:\n", ITvector.size(), mass);

  vcg::Matrix33f IT;
  vcg::Point3f COM;
  vcg::Point3f Trans;
  double link_mass;
  for (size_t i = 0; i < ITvector.size(); i++)
  {
      printf("%s:\n", argv[i+1]);
      IT = ITvector.at(i) * static_cast<float>(mass / volume);
      COM = COMvector.at(i);
      Trans.SetZero();
      for (size_t j = 0; j < joints.size() && j <= i; j++) Trans -= vcg::Point3f(joints.at(j)[0], joints.at(j)[1], joints.at(j)[2]);
      COM += Trans;
      link_mass = mass * VolPerLink.at(i) / volume;
      printf("        <inertial>\n");
      printf("            <mass value=\"%f\" />\n", link_mass);
      printf("            <origin rpy=\"0 0 0\" xyz=\"%014.11f %014.11f %014.11f\" />\n", COM[0], COM[1], COM[2]);
      printf("            <inertia ixx=\"%014.11f\" ixy=\"%014.11f\" ixz=\"%014.11f\"\n",   IT[0][0], IT[0][1], IT[0][2]);
      printf("                                          iyy=\"%014.11f\" iyz=\"%014.11f\"\n", IT[1][1], IT[1][2]);
      printf("                                                               izz=\"%014.11f\" />\n", IT[2][2]);
      printf("        </inertial>\n");
      printf("        <visual>\n");
      printf("            <origin rpy=\"0 0 0\" xyz=\"%014.11f %014.11f %014.11f\" />\n", Trans[0], Trans[1], Trans[2]);
      printf("            <geometry>\n");
      printf("                <mesh filename=\"model://%s\" />\n", argv[i+1]);
      printf("            </geometry>\n");
      printf("        </visual>\n");
  }


  return 0;
}
