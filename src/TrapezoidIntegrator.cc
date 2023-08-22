//---------------------------------------------------------------------------------------------------------
// Module:      TrapezoidIntegrator.c
// Description: Implements a numerical timestep integrator using
//              the Trapezoid integration method.
//
// Includes private functions:
//     Nodes_Advance()
//---------------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"
#include "Node.h"
#include "V3.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

#define DIR_TIMESTEP_ERROR "timestep_error"
#define DIR_NODAL_TIMESTEP "timestep"

#ifdef DEBUG_TIMESTEP_ERROR
//---------------------------------------------------------------------------------------------------------
// Function:    DumpTimestepError
// Description: This is a debug-only funcion which can write
//              to files the per-node positioning error
//              encountered during the timestep integration.
//
//              Positioning error data is written once each
//              timestep with a separate file for each domain.
//              The format of each data line of the files is:
//
//                  <error> <attempted_deltaT> # (<node_id>)
//
//              Files are created with names of the format:
//
//                  timestep_error/NNNNNN/MMMMMM
//
//              where NNNNNN is the cycle number and MMMMMM is
//              the domain ID.
//
// Arguments:
//     deltaT   Duration of the attempted timestep
//---------------------------------------------------------------------------------------------------------

static void DumpTimestepError(Home_t *home, real8 deltaT)
{
        Param_t *param = home->param;

        // Set directory and file names under which to put the
        // data files for this cycle.  Then have task zero
        // create the necessary directories if they don't
        // exist.  Once the directories are in place, all tasks
        // write their own files.

        char  subdir  [256]; snprintf(subdir  , sizeof(subdir  ), "./%s"   , DIR_TIMESTEP_ERROR      );
        char  stepdir [256]; snprintf(stepdir , sizeof(stepdir ), "%s/%06d", subdir , home->cycle+1  );
        char  filename[256]; snprintf(filename, sizeof(filename), "%s/%06d", stepdir, home->myDomain );

        if (home->myDomain == 0)
        {
            (void) mkdir(subdir , S_IRWXU);
            (void) mkdir(stepdir, S_IRWXU);
        }

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        FILE *fd = fopen(filename, "w");

        if (!fd) { Fatal("Open error on timestep error data file %s", filename); }

        for (int i=0; i<home->newNodeKeyPtr; i++)
        {
            Node_t *node = home->nodeKeys[i];

            if (node)
            {
               real8 err = 0.0;

               real8 oldx = node->oldx;
               real8 oldy = node->oldy;
               real8 oldz = node->oldz;

               PBCPOSITION(param, node->x, node->y, node->z, &oldx, &oldy, &oldz);

               err = MAX(err, fabs(node->x - oldx - ((node->vX+node->oldvX)*0.5*deltaT)));
               err = MAX(err, fabs(node->y - oldy - ((node->vY+node->oldvY)*0.5*deltaT)));
               err = MAX(err, fabs(node->z - oldz - ((node->vZ+node->oldvZ)*0.5*deltaT)));

               fprintf(fd, "%14e  %14e  # (%d,%d)\n", err, deltaT, node->myTag.domainID, node->myTag.index);
            }
        }

        fclose(fd);

        return;
}
#endif

#ifdef DEBUG_NODAL_TIMESTEP
//---------------------------------------------------------------------------------------------------------
// Function:    DumpPerNodeTimestep
// Description: This is a debug-only funcion which can write
//              to files the estimated  maximum deltaT
//              (per-node) that could have been taken
//              based on the current velocity and nodal
//              position changes.  This estimate is based
//              on the assumption that each node will move
//              with no more positioning error than the
//              error tolerance specified by <param->rTol>.
//
//              Max deltaT data is written once each
//              timestep with a separate file for each domain.
//              The format of each data line of the files is:
//
//                  <maxDT> <currentDT> # (<node_id>)
//
//              Files are created with names of the format:
//
//                  timestep/NNNNNN/MMMMMM
//
//              where NNNNNN is the cycle number and MMMMMM is
//              the domain ID.
//
// Arguments:
//     newDT   The deltaT selected for this timestep.
//---------------------------------------------------------------------------------------------------------

static void DumpPerNodeTimestep(Home_t *home, double newDT)
{
        real8   signFact;
        real8   errMax, xErr, yErr, zErr;
        real8   oldX, oldY, oldZ;
        real8   newX, newY, newZ;
        real8   oldVX, oldVY, oldVZ;
        real8   newVX, newVY, newVZ;
        real8   xDTMax, yDTMax, zDTMax, nodalMaxDT, nodalMaxErr, minDT;
        Node_t  *node;
        Param_t *param;

        param = home->param;
        minDT = param->maxDT;

        // Set directory and file names under which to put the
        // data files for this cycle.  Then have task zero
        // create the necessary directories if they don't
        // exist.  Once the directories are in place, all tasks
        // write their own files.

        char subdir  [256]; snprintf(subdir  , sizeof(subdir  ), "./%s"   , DIR_NODAL_TIMESTEP     );
        char stepdir [256]; snprintf(stepdir , sizeof(stepdir ), "%s/%06d", subdir, home->cycle+1  );
        char filename[256]; snprintf(filename, sizeof(filename), "%s/%06d", stepdir, home->myDomain);

        if (home->myDomain == 0)
        {
            (void) mkdir(subdir, S_IRWXU);
            (void) mkdir(stepdir, S_IRWXU);
        }

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        FILE *fd = fopen(filename,"w");

        if (!fd) { Fatal("Open error on nodal timestep data file %s", filename); }

        for (int i=0; i < home->newNodeKeyPtr; i++)
        {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;

            // If the node is pinned in all dimensions, it is immobile, so
            // just skip it.

            if (HAS_ALL_OF_CONSTRAINTS(node->constraint, PINNED_NODE)) {
                continue;
            }

            newX = node->x;    newY = node->y;    newZ = node->z;
            oldX = node->oldx; oldY = node->oldy; oldZ = node->oldz;

            PBCPOSITION(param, newX, newY, newZ, &oldX, &oldY, &oldZ);

            newVX = node->vX;    newVY = node->vY;    newVZ = node->vZ;
            oldVX = node->oldvX; oldVY = node->oldvY; oldVZ = node->oldvZ;

            xErr = newX - oldX - ((newVX+oldVX)*0.5*newDT);
            yErr = newY - oldY - ((newVY+oldVY)*0.5*newDT);
            zErr = newZ - oldZ - ((newVZ+oldVZ)*0.5*newDT);

            errMax = param->rTol;

            signFact = ((oldVX+newVX) < 0.0 ? -1.0 : 1.0);
            xDTMax = (2.0 * ((errMax * signFact) + newX - oldX)) / (oldVX + newVX);

            signFact = ((oldVY+newVY) < 0.0 ? -1.0 : 1.0);
            yDTMax = (2.0 * ((errMax * signFact) + newY - oldY)) / (oldVY + newVY);

            signFact = ((oldVZ+newVZ) < 0.0 ? -1.0 : 1.0);
            zDTMax = (2.0 * ((errMax * signFact) + newZ - oldZ)) / (oldVZ + newVZ);

            nodalMaxErr = MAX(MAX(fabs(xErr), fabs(yErr)), fabs(zErr));
            nodalMaxDT  = MIN(MIN(xDTMax, yDTMax), zDTMax);

            if (nodalMaxDT < minDT) minDT = nodalMaxDT;

            fprintf(fd, "%14e  %14e   # (%d,%d)\n", nodalMaxDT,
                    newDT, node->myTag.domainID, node->myTag.index);
        }

#if 0
        printf(" Minimum DeltaT = %e\n", minDT);
#endif

        fclose(fd);

        return;
}
#endif


#ifdef DEBUG_TIMESTEP
//---------------------------------------------------------------------------------------------------------
// DumpNode()
//
// The following is just a debug function used to dump information about the nodes that are
// causing the timestep to be cut.
//---------------------------------------------------------------------------------------------------------

static void DumpNode(Node_t *node)
{
        if (node == (Node_t *)NULL) return;

#if 1
        // print the nodal tag and neighbor tags

        printf("  node(%d,%d) arms %d,  ",
               node->myTag.domainID, node->myTag.index, node->numNbrs);

        for (int i=0; i < node->numNbrs; i++)
            printf(" (%d,%d)", node->nbrTag[i].domainID, node->nbrTag[i].index);

        printf("\n");
#endif

#if 1
        // print the nodal velocity and total node force

        printf("  node(%d,%d)     position = (%e %e %e)\n", node->myTag.domainID, node->myTag.index, node->x, node->y, node->z);
        printf("  node(%d,%d) old position = (%e %e %e)\n", node->myTag.domainID, node->myTag.index, node->oldx, node->oldy, node->oldz);
#endif

#if 1
        // print the nodal velocity and total node force

        printf("  node(%d,%d)    v = (%e %e %e)\n", node->myTag.domainID, node->myTag.index, node->vX, node->vY, node->vZ);
        printf("  node(%d,%d)    f = (%e %e %e)\n", node->myTag.domainID, node->myTag.index, node->fX, node->fY, node->fZ);
#endif

#if 0
        // print the old nodal velocity and total node force

        printf("  node(%d,%d) oldv = (%e %e %e)\n", node->myTag.domainID, node->myTag.index, node->oldvX, node->oldvY, node->oldvZ);
#endif

#if 0
        // print the arm specific forces

        for (int i=0; i < node->numNbrs; i++)
        {
            printf("  node(%d %d) arm[%d]-> (%d %d) f = (%e %e %e)\n",
                   node->myTag.domainID, node->myTag.index, i,
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->armfx[i], node->armfy[i], node->armfz[i]);
        }
#endif

#if 0
        // print the burger's vector for each arm of the node

        for (int i=0; i < node->numNbrs; i++)
        {
            printf("  node(%d %d) arm[%d]-> (%d %d) b = (%f %f %f)\n",
                   node->myTag.domainID, node->myTag.index, i,
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->burgX[i], node->burgY[i], node->burgZ[i]);
        }
#endif

#if 0
        // print the glide plane normal for each arm of the node

        for (int i=0; i < node->numNbrs; i++)
        {
            printf("  node(%d %d) arm[%d]-> (%d %d) n = (%f %f %f)\n",
                   node->myTag.domainID,node->myTag.index, i,
                   node->nbrTag[i].domainID,node->nbrTag[i].index,
                   node->nx[i],node->ny[i],node->nz[i]);
        }
#endif

        return;
}
#endif  // if DEBUG_TIMESTEP


// Limit_Node_Velocity()
//
// It's possible that velocities computed inside the mobility functions can result in non-physical
// node velocities. This routine will limit the computed velocities within this module to the
// shear sound velocity provided via param.
//
// If the shear velocity is greater than zero and the magnitiude of the velocity vector exceeds
// the shear velocity, the velocity is rescaled to a shear velocity vector directed along the
// original vector.
//---------------------------------------------------------------------------------------------------------

#ifdef LIMIT_NODE_VELOCITY
static void Limit_Node_Velocity
(
         real8  & vx  ,   ///< input velocity (x) (b/s) (potentially limited)
         real8  & vy  ,   ///< input velocity (y) (b/s) (potentially limited)
         real8  & vz  ,   ///< input velocity (z) (b/s) (potentially limited)
   const real8    bmag,   ///< magnitude of the burgers vector (m)
   const real8    vmax    ///< shear sound velocity (m/s)
)
{
   if (vmax>0.0)
   {
      real8 vmag = sqrt(vx*vx + vy*vy + vz*vz) * bmag;   // vmag = magnitude of node velocity (m/s)
      real8 vs   = ( (vmag>vmax) ? (vmax/vmag) : 1.0 );  // vs   = scale factor

            vx *= vs;                                    // apply scale factor
            vy *= vs;                                    //
            vz *= vs;                                    //
   }
}
#else
#define Limit_Node_Velocity(a,b,c,d,e)
#endif

//---------------------------------------------------------------------------------------------------------
// Advance_Node()
//
// Will update the node's position given a timestep (dt). This method will average
// the current and previous node velocities for a given timestep.
//
// Additional parameters are provided to adjust for periodic boundary conditions (PBC)...
//
//    cx,cy,cz - the center of the simulation box <xyz>
//    lx,ly,lz - the size   of the simulation box <xyz>
//    sx,sy,sz - PBC adjustment reciprocals       <xyz>
//                 If PBC for a particular axis is active, the PBC reciprocal will be
//                 the reciprocal of the simulation size, otherwise zero.
//                 e.g  sx = ( pbcx && (lx>0.0) ? (1.0/lx) : 0.0 );
//                      sy = ( pbcy && (ly>0.0) ? (1.0/ly) : 0.0 );
//                      sz = ( pbcz && (lz>0.0) ? (1.0/lz) : 0.0 );
//                 By setting  sx,sy,sz to the reciprocals, the rint() correction will be applied.
//---------------------------------------------------------------------------------------------------------

static void Advance_Node
(
         Node_t *node,                               ///< node = pointer to a node
   const real8 dt,                                   ///< dt = next time step
   const real8 cx, const real8 cy, const real8 cz,   ///< simulation center <xyz>
   const real8 lx, const real8 ly, const real8 lz,   ///< simulation size   <xyz>
   const real8 sx, const real8 sy, const real8 sz,   ///< pbc adjustment reciprocals <xyz>
   const real8 bmag,                                 ///< magnitude of the burgers vector (m)
   const real8 vmax,                                 ///< shear sound velocity            (m/s)
   const real8 pmax                                  ///< max allowed position shift (typically maxSeg) (b)
)
{
   if (node)
   {
      // If we don't have a value for the previous velocity, assume
      // previous velocity was same as current velocity.

      if (    (node->oldvX == 0.0)
           && (node->oldvY == 0.0)
           && (node->oldvZ == 0.0) )
      {
         node->oldvX = node->currvX;
         node->oldvY = node->currvY;
         node->oldvZ = node->currvZ;
      }

      // compute the velocity...

      real8 vx = (node->currvX + node->oldvX)/2.0;  // vx = average of current and previous velocity (x)
      real8 vy = (node->currvY + node->oldvY)/2.0;  // vy = average of current and previous velocity (y)
      real8 vz = (node->currvZ + node->oldvZ)/2.0;  // vz = average of current and previous velocity (z)

      Limit_Node_Velocity(vx,vy,vz, bmag,vmax);     // limit velocity to shear sound velocity

      real8 dx = (vx*dt);
      real8 dy = (vy*dt);
      real8 dz = (vz*dt);

#if 0
      real8 dp = sqrt( dx*dx + dy*dy + dz*dz);

      if (dp>pmax)
      {
         printf("%s::%s(%d) WARNING - large node displacement in integration node=(%2d,%4d), dt=%le dpos=%6.1lf pmax=%6.1lf\n",
                __FILE__, __func__, __LINE__,
                node->myTag.domainID, node->myTag.index, dt, dp, pmax );
      }
#endif

      // update the position...

      real8 x = node->oldx + dx;
      real8 y = node->oldy + dy;
      real8 z = node->oldz + dz;

      // adjust for periodic boundary and copy back into node...

      node->x = x - rint((x-cx)*sx)*lx;
      node->y = y - rint((y-cy)*sy)*ly;
      node->z = z - rint((z-cz)*sz)*lz;
   }
}

//---------------------------------------------------------------------------------------------------------
// Advance_Nodes()
//
// Advances an array of nodes or node pointers. Basically - sets up the parallel for-loop that
// drives the nodal advance.
//---------------------------------------------------------------------------------------------------------

static void Advance_Nodes
(
         Node_t **narr,                                 ///< array of node pointers (can be sparse)
   const int      ncnt,                                 ///< number of nodes
   const real8    dt,                                   ///< dt = next time step
   const real8    cx, const real8 cy, const real8 cz,   ///< simulation center <xyz>
   const real8    lx, const real8 ly, const real8 lz,   ///< simulation size   <xyz>
   const real8    sx, const real8 sy, const real8 sz,   ///< pbc adjustment reciprocals <xyz>
   const real8    bmag,                                 ///< magnitude of the burgers vector (m)
   const real8    vmax,                                 ///< shear sound velocity            (m/s)
   const real8    pmax                                  ///< max allowed position shift (typically maxSeg) (b)
)
{
   if (narr && (ncnt>0))
   {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i=0; (i<ncnt); i++)
         Advance_Node(narr[i], dt, cx,cy,cz, lx,ly,lz, sx,sy,sz, bmag,vmax,pmax);
   }
}

//---------------------------------------------------------------------------------------------------------
//  Function:     Advance_Nodes
//  Description:  Advances the position of all nodes (local and ghost) based
//                on the old/current nodal velocities and time deltas.
//
//                This function assumes current nodal positions and velocities
//                have already been preserved in the old* variables by a call
//                to PreserveNodalData().
//
//  Given the following:
//
//      currDT = desired delta time this timestep
//      oldDT  = delta time used in the previous timestep
//      currV  = velocity on entry to timestep integrator
//      oldV   = velocity on entry to timestep integrator on previous step
//      currP  = current nodal position
//      newP   = initial guess at position node will end up in after this timestep.
//
//      The new positions are calculates as:
//
//          newP = currP + 0.5 * (currV + oldV) * currDT;
//---------------------------------------------------------------------------------------------------------

static void Advance_Nodes(Home_t *home, real8 oldDT)
{
   Param_t *param = home->param;

   param->realdt = param->deltaTT;

   const real8 dt    = param->realdt;
   const real8 bmag  = param->burgMag;           // bmag = magnitude of the burgers vector (m)
   const real8 vmax  = param->shearVelocity;     // vmax = shear sound velocity (m/s)
   const real8 pmax  = param->maxSeg;            // pmax = max allowed position shift (b)

   // set simlulation center...

   const real8 cx = (param->maxSideX + param->minSideX)/2.0;
   const real8 cy = (param->maxSideY + param->minSideY)/2.0;
   const real8 cz = (param->maxSideZ + param->minSideZ)/2.0;

   // set simulation size...

   const real8 lx = (param->maxSideX - param->minSideX);
   const real8 ly = (param->maxSideY - param->minSideY);
   const real8 lz = (param->maxSideZ - param->minSideZ);

   // set simulation periodic boundary reciprocals...

   const real8 sx = ( (param->xBoundType == Periodic) && (lx>0.0) ? (1.0/lx) : 0.0 );
   const real8 sy = ( (param->yBoundType == Periodic) && (ly>0.0) ? (1.0/ly) : 0.0 );
   const real8 sz = ( (param->zBoundType == Periodic) && (lz>0.0) ? (1.0/lz) : 0.0 );

   // advance the native and ghost node arrays...

   Advance_Nodes(home->nodeKeys     , home->newNodeKeyPtr , dt, cx,cy,cz, lx,ly,lz, sx,sy,sz, bmag,vmax,pmax);
   Advance_Nodes(home->ghostNodeList, home->ghostNodeCount, dt, cx,cy,cz, lx,ly,lz, sx,sy,sz, bmag,vmax,pmax);
}

//---------------------------------------------------------------------------------------------------------
// Reposition_Node()
//
// Given a pointer to a node, repositions a node for the next iteration of the solver.
//---------------------------------------------------------------------------------------------------------

static void Reposition_Node
(
         Node_t  *node ,  ///< node = pointer to a node
   const real8    dt   ,  ///< dt   = timestep
   const Param_t *param,  ///< param = points to the current param data structure
   const real8    bmag ,  ///< magnitude of the burgers vector (m)
   const real8    vmax    ///< shear sound velocity            (m/s)
)
{
   if (node)
   {
      real8 x  = node->x;                        // x  = current  position (x)
      real8 y  = node->y;                        // y  = current  position (y)
      real8 z  = node->z;                        // z  = current  position (z)

      real8 px = node->oldx;                     // px = previous position (x)
      real8 py = node->oldy;                     // py = previous position (y)
      real8 pz = node->oldz;                     // pz = previous position (z)

      real8 vx = (node->vX + node->currvX)/2.0;  // vx = average  velocity (x)
      real8 vy = (node->vY + node->currvY)/2.0;  // vy = average  velocity (y)
      real8 vz = (node->vZ + node->currvZ)/2.0;  // vz = average  velocity (z)

      Limit_Node_Velocity(vx,vy,vz, bmag,vmax);  // limit velocity to shear sound velocity

      real8 dx = (vx*dt);                        // dx = change in position (x)
      real8 dy = (vy*dt);                        // dy = change in position (y)
      real8 dz = (vz*dt);                        // dz = change in position (z)

#if 0
      // generate diagnostic on large nodal moves...

      real8 dp = sqrt(dx*dx + dy*dy + dz*dz);    // dp = magnitide of change in position

      if (dp > param->maxSeg )
      {
         Tag_t tag   = node->myTag;

         real8 v [3] = { vx          , vy          , vz           };  V3_SCALE(v ,v ,bmag);   // v  = average  velocity (m/s)
         real8 v0[3] = { node->vX    , node->vY    , node->vZ     };  V3_SCALE(v0,v0,bmag);   // v0 = previous velocity (m/s)
         real8 v1[3] = { node->currvX, node->currvY, node->currvZ };  V3_SCALE(v1,v1,bmag);   // v1 = current  velocity (m/s)

         real8 vm    = V3_MAG(v );  // vm  = magnitude of average  velocity (m/s)
         real8 v0m   = V3_MAG(v0);  // v0m = magnitude of previous velocity (m/s)
         real8 v1m   = V3_MAG(v1);  // v1m = magnitude of current  velocity (m/s)

         printf("%s::%s(%d) WARNING - large node displacement in integration node=(%2d,%4d), dt=%le dp=%6.1lf vm=%6.1lf vshear=%6.1lf vXYZ=%6.1lf  [%6.1lf %6.1lf %6.1lf] currvXYZ=%6.1lf [%6.1lf %6.1lf %6.1lf]\n",
                __FILE__, __func__, __LINE__,
                tag.domainID, tag.index,
                dt, dp,
                vm, vmax,
                v0m, v0[0],v0[1],v0[2],
                v1m, v1[0],v1[1],v1[2] );
      }
#endif

      PBCPOSITION(param, x,y,z, &px, &py, &pz);

      x = px + dx;
      y = py + dy;
      z = pz + dz;

      FoldBox(param, &x, &y, &z);

      node->x = x;
      node->y = y;
      node->z = z;
   }
}

//---------------------------------------------------------------------------------------------------------
// Reposition_Nodes()
//
// Executes a parallel reposition of an array of node pointers.
// Returns the magnitude of the distance of the node with the largest change in position
//---------------------------------------------------------------------------------------------------------

static void Reposition_Nodes
(
         Node_t **narr ,  ///< array of node pointers (can be sparse)
   const int      ncnt ,  ///< number of nodes
   const real8    dt   ,  ///< next time step
   const Param_t *param   ///< points to the current param data structure
)
{
   if (narr && (ncnt>0))
   {
      const real8 bmag  = param->burgMag;           // bmag = magnitude of the burgers vector (m)
      const real8 vmax  = param->shearVelocity;     // vmax = shear sound velocity (m/s)

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i=0; (i<ncnt); i++)
      {
         Reposition_Node(narr[i],dt,param,bmag,vmax);
      }
   }
}

//---------------------------------------------------------------------------------------------------------
// Preserve_Node()
//
// Will preserve various node data elements specified via the items argument.
//
// The items argument is a bitwise 'or' of the following constants...
//    NODE_POSITION : saves the current node position into the old position
//    NODE_CURR_VEL : saves the current velocity without wiping out the copy of the previous velocity
//    NODE_OLD_VEL  : saves the old velocity (normally done at the end of the timestep integrator
//---------------------------------------------------------------------------------------------------------

static void Preserve_Node (Node_t *node, const int items)
{
   if (node)
   {
      if (items & NODE_POSITION)
      {
         node->oldx   = node->x;
         node->oldy   = node->y;
         node->oldz   = node->z;
      }

      if (items & NODE_CURR_VEL)
      {
         node->currvX = node->vX;
         node->currvY = node->vY;
         node->currvZ = node->vZ;
      }

      if (items & NODE_OLD_VEL)
      {
         node->oldvX  = node->currvX;
         node->oldvY  = node->currvY;
         node->oldvZ  = node->currvZ;
      }
   }
}

//---------------------------------------------------------------------------------------------------------
// Preserve_Nodes()
//
// Preserves the elements of an array of node pointers.
//---------------------------------------------------------------------------------------------------------

static void Preserve_Nodes
(
         Node_t **narr ,  ///< array of node pointers (can be sparse)
   const int      ncnt ,  ///< number of nodes
   const int      items   ///< items to preserve (see above)
)
{
   if (narr && (ncnt>0))
   {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i=0; (i<ncnt); i++)
         Preserve_Node(narr[i],items);
   }
}

static real8 Max(real8 a, real8 b, real8 c)
{
   return ( (a>b) ? ( (a>c) ? a : c ) : ( (b>c) ? b : c ) );
}

//---------------------------------------------------------------------------------------------------------
// Node_Perror()
//
// Computes the maximum position error of the current node and a given timestep.
//
// The returned error is the maximum of the individual vector component errors.
//---------------------------------------------------------------------------------------------------------

static real8 Node_Perror
(
   const Node_t  *node ,    ///< points to node
   const Param_t *param,    ///< points to param structure
   const real8    dt   ,    ///< current time step
   const real8    bmag ,    ///< bmag = magnitude of the burgers vector (m)
   const real8    vmax      ///< vmax = shear sound velocity (m/s)
)
{
    real8 perr=0.0;

    if (node)
    {
        real8 x    = node->x;                      // x  = current  position (x)
        real8 y    = node->y;                      // y  = current  position (y)
        real8 z    = node->z;                      // z  = current  position (z)

        real8 px   = node->oldx;                   // px = previous position (x)
        real8 py   = node->oldy;                   // py = previous position (y)
        real8 pz   = node->oldz;                   // pz = previous position (z)

        PBCPOSITION(param, x,y,z, &px, &py, &pz);  // adjust previous position to nearest PBC image to current position

        real8 vx   = (node->vX+node->currvX)/2.0;  // vx = average of current and previous velocity (x)
        real8 vy   = (node->vY+node->currvY)/2.0;  // vy = average of current and previous velocity (y)
        real8 vz   = (node->vZ+node->currvZ)/2.0;  // vz = average of current and previous velocity (z)

        Limit_Node_Velocity(vx,vy,vz, bmag,vmax);  // limit velocity to shear sound velocity

        real8 dx = (vx*dt);                        // dx = change in position (x)
        real8 dy = (vy*dt);                        // dy = change in position (y)
        real8 dz = (vz*dt);                        // dz = change in position (z)

        perr = Max( fabs(x-px-dx),                 // perr = max of all position errors (xyz)
                    fabs(y-py-dy),                 //
                    fabs(z-pz-dz) );               //
    }

    return(perr);
}

//---------------------------------------------------------------------------------------------------------
// Node_dPos()
//
// Returns the magnitude of the current position change given velocities and time step.
//---------------------------------------------------------------------------------------------------------

static real8 Node_dPos
(
   const Node_t  *node  ,   ///< points to node
   const real8    dt    ,   ///< current time step
   const real8    bmag  ,   ///< magnitude of the burgers vector (m)
   const real8    vmax      ///< shear sound velocity            (m/s)
)
{
    real8 dp=0.0;

    if (node)
    {
        real8 vx = (node->vX+node->currvX)/2.0;    // vx = average of previous and current velocity (x)
        real8 vy = (node->vY+node->currvY)/2.0;    // vy = average of previous and current velocity (y)
        real8 vz = (node->vZ+node->currvZ)/2.0;    // vz = average of previous and current velocity (z)

        Limit_Node_Velocity(vx,vy,vz, bmag,vmax);  // limit velocity to shear sound velocity

        real8 dx = (vx*dt);                        // dx = change in position (x)
        real8 dy = (vy*dt);                        // dy = change in position (y)
        real8 dz = (vz*dt);                        // dz = change in position (z)

              dp = sqrt( dx*dx + dy*dy + dz*dz );  // dp = magnitude of position change
    }

    return(dp);
}

//---------------------------------------------------------------------------------------------------------
//  Function:    TrapezoidIntegrator
//  Description: Implements a numerical timestep integrator using
//               the Trapezoid integration method.
//
//  Note: This function assumes that the nodal force/velocity data is accurate
//  for the current positions of the nodes on entry to the routine.
//---------------------------------------------------------------------------------------------------------

void TrapezoidIntegrator(Home_t *home)
{

#ifdef DEBUG_TIMESTEP_ERROR
        int     dumpErrorData=1;
#endif

        Node_t  *maxNode = 0;                 // maxNode = points to node with largest position error
        Node_t  *posNode = 0;                 // posNode = points to node with largest change in position

        Param_t *param   = home->param;

        real8    oldDT   = param->deltaTT;
        real8    newDT   = MIN(param->maxDT, param->nextDT);

        if (newDT <= 0.0) newDT = param->maxDT;

        param->deltaTT = newDT;

        // Preserve certain items of nodal data for later use.

        Preserve_Nodes (home->nodeKeys     , home->newNodeKeyPtr , NODE_POSITION | NODE_CURR_VEL );
        Preserve_Nodes (home->ghostNodeList, home->ghostNodeCount, NODE_POSITION | NODE_CURR_VEL );

        // Loop until we converge on a time step.  First step is to
        // use the current positions and velocities andreposition the
        // nodes to where they would be after the suggested delta time.

        int   convergent      = 0;
        int   incrDelta       = 1;
        int   maxIterations   = param->trapezoidMaxIterations;

        real8 errMax          = 0.0;
        real8 posMax          = 0.0;
        real8 globalErrMax    = 0.0;
        real8 globalPosMax    = 0.0;
        int   globalIterError = 0;

        while (!convergent)
        {
            // Advance all nodes from their previous positions to new positions
            // based on their current velocities and the delta T being tested.
            // This includes all ghost nodes, so nodal velocities must previously
            // have been distributed to neighboring domains.

            Advance_Nodes(home, oldDT);

            int mobIterError=0;

            for (int iter=0; iter < maxIterations; iter++)
            {
                // Recalculate nodal force and velocity at the new positions

                int  doAll=1;
                NodeForce(home, FULL);
                mobIterError = CalcNodeVelocities(home, 0, doAll);
                CommSendVelocity(home);

                // Note - there's no need to compute position errors if the
                // mobility function didn't converge...

                errMax=0.0;
                posMax=0.0;

                if (!mobIterError)
                {
                    // compute position errors...

                    for (int i=0; i<home->newNodeKeyPtr; i++)
                    {
                        Node_t *node = home->nodeKeys[i];

                        if (node)
                        {
                           const real8 bmag = param->burgMag;                      // bmag = magnitude of the burgers vector (m)
                           const real8 vmax = param->shearVelocity;                // vmax = shear sound velocity (m/s)

                           real8 perr = Node_Perror(node,param,newDT, bmag,vmax);  // perr = position error

                           maxNode = ( (perr>errMax) ? node : maxNode );
                           errMax  = ( (perr>errMax) ? perr : errMax  );

                                 real8 dpos = Node_dPos(node,newDT, bmag,vmax);    // dpos = change in position

                           posNode = ( (dpos>posMax) ? node : posNode );
                           posMax  = ( (dpos>posMax) ? dpos : posMax  );
                        }
                    }
                }

                // Determine global error across ALL domains...

#ifdef PARALLEL
                real8 localVals [3] = { 0.0 };
                real8 globalVals[3] = { 0.0 };

                localVals[0]    = errMax;
                localVals[1]    = posMax;
                localVals[2]    = (real8) mobIterError;

                MPI_Allreduce(localVals, globalVals, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

                globalErrMax    = globalVals[0];
                globalPosMax    = globalVals[1];
                globalIterError = globalVals[2];
#else
                globalErrMax    = errMax;
                globalPosMax    = posMax;
                globalIterError = mobIterError;
#endif

                // If any domain encountered an error iterating inside
                // the mobility function, just go right to cutting the
                // timestep and trying again.

                if (globalIterError) {
                    iter = maxIterations;
                    continue;
                }

                // If the error is within the tolerance AND the maximum node position change is less than maxSeg,
                // we've reached convergence so we can accept this deltaT. Otherwise reposition the nodes and try again.

                if (    (globalErrMax < param->rTol  )
                     && (globalPosMax < param->maxSeg) )
                {
                    convergent = 1;  // (converged)
                    break;
                }
                else
                {
                    convergent = 0; // (not converged)
                    incrDelta  = 0;

                    if (iter == maxIterations-1) { continue; }

                    Reposition_Nodes(home->nodeKeys     , home->newNodeKeyPtr , newDT, param);  // reposition native nodes
                    Reposition_Nodes(home->ghostNodeList, home->ghostNodeCount, newDT, param);  // reposition ghost  nodes
                }
            }  // for (iter = 0; ...)

#ifdef DEBUG_TIMESTEP_ERROR
            if (dumpErrorData)
            {
                DumpTimestepError(home, newDT);
                dumpErrorData = 0;
            }
#endif

            // If we haven't converged, cut the timestep and try again...

            if (!convergent)
            {
                newDT *= param->dtDecrementFact;
                param->deltaTT = newDT;

                if ((newDT < 1.0e-20) && (home->myDomain == 0)) {
                    Fatal("TrapezoidIntegrator(): Timestep has dropped below minimal threshold to %e.  Aborting!", newDT);
                }

#ifdef DEBUG_TIMESTEP
                if ((home->myDomain == 0) && (globalIterError))
                {
                    printf(" +++ Cut timestep to %e for mobility non-convergence\n", newDT);
                }

                // If this is this domain with the node causing timestep to drop
                // dump some info on the node limiting the timestep.

                if ((globalErrMax == errMax) && (globalIterError == 0) && maxNode )
                {
                    printf("  Cut timestep for (%d,%d):  errMax %e, newDT %e\n", maxNode->myTag.domainID, maxNode->myTag.index, errMax, newDT);
                    DumpNode(maxNode);
                }
#endif
            }  // end if (!convergent)

        }  // while (!convergent)

#ifdef DEBUG_NODAL_TIMESTEP
        DumpPerNodeTimestep(home, newDT);
#endif

        // Automatically increment timestep if convergence was reached
        // on the very first iteration of the above loop.
        //
        // If variable timestep adjustments are enabled, calculate an
        // adjustment factor based on the maximum allowed timestep increment
        // and the maximum error found above.  If variable timestep
        // adjustments are not enabled, adjust the timestep by the
        // maximum permitted factor.

        param->deltaTT   = newDT;
        param->realdt    = newDT;
        param->timeStart = param->timeNow;

        if (incrDelta)
        {
            if (param->dtVariableAdjustment)
            {
                real8 tmp1   = pow(param->dtIncrementFact, param->dtExponent);
                real8 tmp2   = globalErrMax/param->rTol;
                real8 tmp3   = 1.0 / param->dtExponent;
                real8 tmp4   = pow(1.0/(1.0+(tmp1-1.0)*tmp2), tmp3);
                real8 factor = param->dtIncrementFact * tmp4;

                param->nextDT = MIN(param->maxDT, newDT*factor);
            }
            else
            {
                param->nextDT = MIN(param->maxDT, newDT*param->dtIncrementFact);
            }
        }
        else
        {
            param->nextDT = newDT;
        }

        // Copy the nodal velocities that existed on entry to the timestep integrator.

        Preserve_Nodes (home->nodeKeys     , home->newNodeKeyPtr , NODE_OLD_VEL );
        Preserve_Nodes (home->ghostNodeList, home->ghostNodeCount, NODE_OLD_VEL );

#ifdef _ARLFEM
        // If we're using the FEM code for simulating free surfaces, we
        // have to handle any nodes/segments that have moved onto or past
        // the simulation surfaces.
        //
        // FIX ME!  Do we need to again recalculate forces/velocities for
        //         the nodes that get modified in FixSurfaceTopology(), or
        //         can we live with the results until the next cycle???

        FixSurfaceTopology(home);
#endif

}
