#!/usr/bin/env pvpython
# or #!/usr/bin/env vtkpython
# e.g. /Programs/Scientific-Visualization/Paraview/ParaView-5.0.0.app/Contents/bin/pvpython
from math import cos, sin, acos, pi, atan2, sqrt
import numpy as np
from numpy import linalg as LA
import vtk, sys, argparse, os

def errexit(msg, code=1):
   print msg
   sys.exit(code)

parser = argparse.ArgumentParser(description='Create a VTK file from an inclusions file')
parser.add_argument('--infile', '-i', default="inclusions.dat", help="name of inclusion file")
parser.add_argument('--outfile', '-o', default="ellipsoids.vtk", help="name of output VTK file")
parser.add_argument('--spherical', '-S', action="store_true", help="Instead of XYZ rotations in fields 8-10, use spherical coordinates (theta and phi) in fields 8 and 9")
parser.add_argument('--help-file-format', "-H", action="store_true", help="Print the accepted file format for the selected input file type (spherical or XYZ rotation)")
args = parser.parse_args()

# =====================================================================
def fileFormatHelp():
   if args.spherical:
      print """# Data for an inclusion consist of the following items separated
# by white space
#
#     Inclusion ID:  Unique integer identifying the inclusion regardless
#                    of its position in the simulation or the domain
#                    encompassing it.  (Ignored)
#     Position:      3 values specifying the X, Y and Z coordinates of
#                    the center of the inclusion
#     Scale:         3 values specifying the X, Y and Z axis lengths before
#                    rotation.
#     theta, phi:    spherical coordinates, theta in XY plane, phi from Z axis
#     Strain field:  (ignored) Strain field for the particle.  Strain field is a
#                    symetric matrix so only six components are specified
#                    in the order: S[0][0], S[1][1], S[2][2], S[1][2],
#                    S[0][2], S[0][1]
#
1   0 0 0    50 5 5  0. 0.  0.1 0.1 0.1 0. 0. 0.
2   0 0 0    50 5 5  90. 90. 0.1 0.1 0.1 0. 0. 0.
3   0 0 0    50 5 5  0. 90. 0.1 0.1 0.1 0. 0. 0."""
   else:
      print """# Data for an inclusion consist of the following items separated
# by white space
#
#     Inclusion ID:  Unique integer identifying the inclusion regardless
#                    of its position in the simulation or the domain
#                    encompassing it.  (Ignored)
#     Position:      3 values specifying the X, Y and Z coordinates of
#                    the center of the inclusion
#     Scale:         3 values specifying the X, Y and Z axis lengths before
#                    rotation.
#     Rotation:      X,Y and Z axis rotations in degrees
#     Strain field:  (ignored) Strain field for the particle.  Strain field is a
#                    symetric matrix so only six components are specified
#                    in the order: S[0][0], S[1][1], S[2][2], S[1][2],
#                    S[0][2], S[0][1]
#
1   0 0 0    50 5 5  0. 0. 0. 0.1 0.1 0.1 0. 0. 0.
1   0 0 0    50 5 5  90 0. 0.  0.1 0.1 0.1 0. 0. 0.
1   0 0 0    50 5 5  0. 90. 0. 0.1 0.1 0.1 0. 0. 0."""
   return

# =====================================================================

if args.help_file_format:
   fileFormatHelp()
   sys.exit(0)

# =====================================================================

if not os.path.exists(args.infile):
   errexit ("File %s does not exist!"%args.infile)

appended = vtk.vtkAppendPolyData()
appended.UserManagedInputsOn()


lines = []
for line in open(args.infile, "r").readlines():
   line = line.strip()
   if line[0] == '#':
      continue
   lines.append(line)

appended.SetNumberOfInputs(len(lines))

idx = 0
for line in lines:
   vals = map(float, line.split())
   ellipsoid = vtk.vtkParametricEllipsoid()

   # ID
   ID = int(vals[0])
   
   # position
   [x,y,z] = [vals[1], vals[2], vals[3]]
   
   # semi principal axis 
   [Sx,Sy,Sz] = [vals[4], vals[5], vals[6]]
   ellipsoid.SetXRadius(Sx)
   ellipsoid.SetYRadius(Sy)
   ellipsoid.SetZRadius(Sz)

   parametricFunctionSource = vtk.vtkParametricFunctionSource()
   parametricFunctionSource.SetParametricFunction(ellipsoid);
   transform = vtk.vtkTransform()
   transformFilter = vtk.vtkTransformFilter()
   transformFilter.SetTransform(transform)
   transformFilter.SetInputConnection(parametricFunctionSource.GetOutputPort())
   transform.Identity()
   transform.Translate(x,y,z)

   # define rotation matrix
   RX = [vals[7], vals[8], vals[9]]
   RX = RX/LA.norm(RX)
   RY = [vals[10], vals[11], vals[12]]
   RY = RY/LA.norm(RY)
   RZ = [vals[8]*vals[12]-vals[9]*vals[11], vals[9]*vals[10]-vals[7]*vals[12], vals[7]*vals[11]-vals[8]*vals[10]];
   RZ = RZ/LA.norm(RZ)

   
   matrix = [RX[0], RY[0], RZ[0], x, 
             RX[1], RY[1], RZ[1], y, 
             RX[2], RY[2], RZ[2], z, 
             0.0, 0.0, 0.0, 1.0]
   transform.SetMatrix(matrix)

   # Apply transformation
   transformFilter.Update()
   sphereoid = transformFilter.GetOutput()
   array = vtk.vtkIntArray()
   array.SetNumberOfComponents(1)
   array.SetName("ellipsoid ID")
   for x in range(0,sphereoid.GetNumberOfPoints()):
      array.InsertNextTuple([ID])
   sphereoid.GetPointData().AddArray(array)
   appended.SetInputDataByNumber(idx, sphereoid)

   idx = idx + 1

writer = vtk.vtkDataSetWriter()
writer.SetFileName(args.outfile)
writer.SetInputConnection(appended.GetOutputPort())
writer.Write()

print "Wrote file", args.outfile

sys.exit(0)


