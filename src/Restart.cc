
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "PBC.h"
#include "Tag.h"
#include "Node.h"
#include "Domain_t.h"
#include "Restart.h"
#include "V3.h"
#include "PBC.h"

//------------------------------------------------------------------------------------------------------------

#define VDELETE(a) if (a) { delete [] a; a=0; }

//---------------------------------------------------------------------------------------------------------

Paradis_Restart::Paradis_Restart
(
   const char *path    ///< path to restart file
)
{
   rs_path          = 0;

   rs_version       = 0;
   rs_nsegs         = 0;

   rs_min_coords[0] = 0.0;
   rs_min_coords[1] = 0.0;
   rs_min_coords[2] = 0.0;

   rs_max_coords[0] = 0.0;
   rs_max_coords[1] = 0.0;
   rs_max_coords[2] = 0.0;

   rs_decomp_type   = 0;
   rs_decomp_geo[0] = 0;
   rs_decomp_geo[1] = 0;
   rs_decomp_geo[2] = 0;

   rs_domains       = 0;
   rs_domain_cnt    = 0;

   rs_nodes         = 0;
   rs_node_cnt      = 0;

   rs_bv_cnt        = 0;
   rs_bv_tbl        = 0;

   Load(path);
}

//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::Recycle (void)
{
   VDELETE(rs_path   );
   VDELETE(rs_domains);
   VDELETE(rs_nodes  );
   VDELETE(rs_bv_tbl );

   rs_domain_cnt = 0;
   rs_node_cnt   = 0;

   rs_bv_cnt     = 0;
   rs_bv_tbl     = 0;
}

//---------------------------------------------------------------------------------------------------------

Paradis_Restart::~Paradis_Restart() { Recycle(); }

//---------------------------------------------------------------------------------------------------------
// Domain sorting...
//---------------------------------------------------------------------------------------------------------

static int dcmp (const void *p0, const void *p1)
{
   Domain_t *d0 = *((Domain_t **) p0);
   Domain_t *d1 = *((Domain_t **) p1);

   if (d0->dndx < d1->dndx ) return(-1);
   if (d0->dndx > d1->dndx ) return(+1);

   return(0);
}

void Paradis_Restart::Sort_Domains (void)
{
   if (rs_domains && (rs_domain_cnt>0))
   {
      Domain_t **darr = new Domain_t *[rs_domain_cnt];

      for (int i=0; (i<rs_domain_cnt); i++)
         darr[i] = rs_domains+i;

      qsort(darr, rs_domain_cnt, sizeof(Domain_t *), dcmp);

      Domain_t *tmp = new Domain_t[rs_domain_cnt];

      for (int i=0; (i<rs_domain_cnt); i++)
         tmp[i] = *darr[i];

      VDELETE(darr);
      VDELETE(rs_domains);

      rs_domains = tmp;
   }
}

//---------------------------------------------------------------------------------------------------------
// Node sorting...
//---------------------------------------------------------------------------------------------------------

static int ncmp (const void *p0, const void *p1)
{
   Node_t *n0 = *((Node_t **) p0);
   Node_t *n1 = *((Node_t **) p1);

   if (n0->myTag.domainID < n1->myTag.domainID) return(-1);   // (n0<n1) (compares domains)
   if (n0->myTag.domainID > n1->myTag.domainID) return(+1);   // (n0>n1) (compares domains)
   if (n0->myTag.index    < n1->myTag.index   ) return(-1);   // (n0<n1) (compares indices)
   if (n0->myTag.index    > n1->myTag.index   ) return(+1);   // (n0>n1) (compares indices)

   return(0);
}

void Paradis_Restart::Sort_Nodes (void)
{
   if (rs_nodes && (rs_node_cnt>0))
   {
      Node_t **narr = new Node_t *[rs_node_cnt];

      for (int i=0; (i<rs_node_cnt); i++)
         narr[i] = rs_nodes+i;

      qsort(narr, rs_node_cnt, sizeof(Node_t *), ncmp);

      Node_t *tmp = new Node_t[rs_node_cnt];

      for (int i=0; (i<rs_node_cnt); i++)
         tmp[i] = *narr[i];

      VDELETE(narr    );
      VDELETE(rs_nodes);

      rs_nodes = tmp;
   }
}

// Find_Node()
//
// Returns a pointer to the node with a given tag. The search is binary and assumes the current
// node list is sorted by tag.
//---------------------------------------------------------------------------------------------------------

Node_t *Paradis_Restart::Find_Node (const Tag_t & tag)
{
   Node_t *node=0;

   if (rs_nodes && (rs_node_cnt>0))
   {
      int i0=0;
      int i1=(rs_node_cnt-1);

      while (i0<=i1)
      {
         int i = (i0+i1)/2;

              node = rs_nodes+i;
         int  cmp  = ( (node->myTag < tag) ? -1 :
                       (node->myTag > tag) ? +1 : 0 );

         if (cmp==0) { return(node); }
         if (cmp< 0) { i0=i+1; }
         else        { i1=i-1; }
      }
   }

   return(node);
}

//---------------------------------------------------------------------------------------------------------
// local file and parsing utilities...
//---------------------------------------------------------------------------------------------------------

static size_t file_size (FILE *fd)
{
   size_t bytes=0;  // file size (in bytes)
   size_t fpos =0;  // saved file position

   if (fd)
   {
      fpos  = ftell(fd);                // save current file position
              fseek(fd,0,SEEK_END);     // seek to the end of the file
      bytes = ftell(fd);                // bytes = file position at end of file
              fseek(fd,fpos,SEEK_SET);  // restore file position
   }

   return(bytes);
}

static char *file_load (const char *path, size_t & bytes)
{
   char *p=0; bytes=0;

   FILE *fd = ( path ? fopen(path,"r") : 0 );

   if (fd)
   {
      bytes = file_size(fd);
      p     = ( (bytes>0) ? new char[bytes] : 0 );

      if (p && (bytes>0))
         fread(p,bytes,1,fd);

      fclose(fd);
   }

   return(p);
}

static char *str_getln(char *p, const char *pe)
{
   if (p && pe && (p<pe))
   {
      while( (p<pe) && (*p!='\r') && (*p!='\n') ) p++;

      if ( (p<pe) && (*p=='\r') ) p++;
      if ( (p<pe) && (*p=='\n') ) p++;
   }

   if (p && (p<pe) && (*p=='#') ) p=str_getln(p,pe);

   return( p && (p<pe) ? p : 0 );
}

static char *str_next(char *p, const char *pe)
{
   if (p && pe)
      while( (p<pe) && isspace(*p) ) p++;

   return( p && (p<pe) ? p : 0 );
}

static int str_match(char *pa, const char *pb)
{
   return ( (pa && pb && (strncmp(pa,pb,strlen(pb))==0)) ? 1 : 0 );
}

//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::Load
(
   const char *path    ///< path to restart file
)
{
   size_t  bytes = 0;
   char   *p     = file_load(path,bytes);
   char   *pe    = ( p ? (p+bytes) : 0 );

   Recycle();

   rs_path = ( path ? strdup(path) : 0 );

   Parse  (p,pe);

   Sort_Domains();
   Sort_Nodes  ();
   Index_Nodes ();

   VDELETE(p);
}

//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::Save (const char *path)
{
   FILE *fd = ( path ? fopen(path,"w") : 0 );

   if (fd) { Print(fd); fclose(fd); }
}

//---------------------------------------------------------------------------------------------------------

static char *Parse_RS_Decomp(char *p, char *pe, Domain_t *doms, int nx, int ny, int nz)
{
   real8   *xb = new real8   [nx+1];
   real8  **yb = new real8  *[nx  ];
   real8 ***zb = new real8 **[nx  ];

   for (int i=0; (i<nx); i++)
   {
      yb[i] = new real8  [ny+1];
      zb[i] = new real8 *[ny  ];

      for (int j=0; (j<ny); j++)
         zb[i][j] = new real8[nz+1];
   }

   int i=0;
   for (i=0; (i<nx); i++)
   {
      p = str_getln(p,pe);
      sscanf(p, "%lf", &xb[i]);

      int j=0;
      for (j=0; (j<ny); j++)
      {
         p = str_getln(p,pe);
         sscanf(p, "%lf", &yb[i][j]);

         int k=0;
         for (k=0; (k<nz); k++)
         {
            p = str_getln(p,pe);
            sscanf(p, "%lf", &zb[i][j][k]);
         }

         p = str_getln(p,pe);
         sscanf(p, "%lf", &zb[i][j][k]);
      }

      p = str_getln(p,pe);
      sscanf(p, "%lf", &yb[i][j]);
   }
   p = str_getln(p,pe);
   sscanf(p, "%lf", &xb[i]);

   int dndx=0;
   for (int i=0; (i<nx); i++)
   for (int j=0; (j<ny); j++)
   for (int k=0; (k<nz); k++, dndx++)
   {
      doms[dndx] = Domain_t( dndx,
                             xb[i]      , xb[i+1],
                             yb[i][j]   , yb[i  ][j+1],
                             zb[i][j][k], zb[i  ][j  ][k+1] );
   }

   for (int i=0; (i<nx); i++)
   {
      for (int j=0; (j<ny); j++)
         delete [] zb[i][j];

      delete [] yb[i];
      delete [] zb[i];
   }

   VDELETE(xb);
   VDELETE(yb);
   VDELETE(zb);

   return(p);
}

static char *Parse_RB_Decomp(char *p, char *pe, Domain_t *doms, int nx, int ny, int nz)
{
   Domain_t *d=doms;

   for (int i=0; (p && (i<(nx*ny*nz))); i++, d++)
   {
      p = str_getln(p,pe);
      sscanf(p,"%d %lf %lf %lf %lf %lf %lf", &d->dndx, &d->x0, &d->y0, &d->z0, &d->x1, &d->y1, &d->z1 );
   }

   return(p);
}

char *Paradis_Restart::Parse_Domains(char *p, char *pe)
{
   VDELETE(rs_domains); rs_domains=0;

   int nx    = rs_decomp_geo[0];
   int ny    = rs_decomp_geo[1];
   int nz    = rs_decomp_geo[2];
   int ndoms = (nx*ny*nz);

   rs_domains    = new Domain_t[ndoms];
   rs_domain_cnt = ndoms;

   if (rs_decomp_type==1) { p = Parse_RS_Decomp(p,pe,rs_domains, nx,ny,nz); }
   if (rs_decomp_type==2) { p = Parse_RB_Decomp(p,pe,rs_domains, nx,ny,nz); }

   return(p);
}

char *Paradis_Restart::Parse_Nodes(char *p, char *pe)
{
   p = str_getln(p,pe);

   VDELETE(rs_nodes);

   rs_nodes    = 0;
   rs_node_cnt = 0;

   int max=0;

   while (p && pe && (p<pe))
   {
      Node_t node;

      p = node.Parse(p,pe);

      rs_nodes = Node_Append(rs_nodes,rs_node_cnt,max,node);
   }

   return(p);
}

//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::Parse
(
   char *p ,   ///< points to start of restart file buffer
   char *pe    ///< points to end   of restart file buffer
)
{
   while (p && pe && (p<pe))
   {
      p = str_next(p,pe);

           if (str_match(p,"dataFileVersion"    )) { sscanf(p,"dataFileVersion = %d"           , &rs_version); }
      else if (str_match(p,"numFileSegments"    )) { sscanf(p,"numFileSegments = %d"           , &rs_nsegs  ); }
      else if (str_match(p,"minCoordinates"     )) { sscanf(p,"minCoordinates = [ %lf %lf %lf" , rs_min_coords, rs_min_coords+1, rs_min_coords+2 ); }
      else if (str_match(p,"maxCoordinates"     )) { sscanf(p,"maxCoordinates = [ %lf %lf %lf" , rs_max_coords, rs_max_coords+1, rs_max_coords+2 ); }
      else if (str_match(p,"dataDecompType"     )) { sscanf(p,"dataDecompType = %d"            , &rs_decomp_type ); }
      else if (str_match(p,"dataDecompGeometry" )) { sscanf(p,"dataDecompGeometry = [ %d %d %d", rs_decomp_geo, rs_decomp_geo+1, rs_decomp_geo+2 ); }

      else if (str_match(p,"domainDecomposition")) { p=Parse_Domains(p,pe); }
      else if (str_match(p,"nodalData"          )) { p=Parse_Nodes  (p,pe); }

      p = str_getln(p,pe);
   }
}


//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::Print (FILE *fd)
{
   if (fd)
   {
      fprintf(fd,"dataFileVersion    = %d\n"                         , rs_version );
      fprintf(fd,"numFileSegments    = %d\n"                         , rs_nsegs   );
      fprintf(fd,"minCoordinates     = [ %12.4lf %12.4lf %12.4lf ]\n", rs_min_coords[0], rs_min_coords[1], rs_min_coords[2] );
      fprintf(fd,"maxCoordinates     = [ %12.4lf %12.4lf %12.4lf ]\n", rs_max_coords[0], rs_max_coords[1], rs_max_coords[2] );
      fprintf(fd,"nodeCount          = %d\n"                         , rs_node_cnt    );
      fprintf(fd,"dataDecompType     = %d\n"                         , 2 );
      fprintf(fd,"dataDecompGeometry = [ %d %d %d ]\n"               , rs_decomp_geo[0], rs_decomp_geo[1], rs_decomp_geo[2] );
      fprintf(fd,"\n");

      if (rs_domains && (rs_domain_cnt>0))
      {
         fprintf(fd,"domainDecomposition = \n");
         fprintf(fd,"# Dom_ID  Minimum XYZ bounds   Maximum XYZ bounds\n");

         Domain_t *dom = rs_domains;

         for (int i=0; (i<rs_domain_cnt); i++, dom++)
            fprintf(fd," %3d %12.4lf %12.4lf %12.4lf %12.4lf %12.4lf %12.4lf\n",
                    dom->dndx, dom->x0, dom->y0, dom->z0, dom->x1, dom->y1, dom->z1 );

         fprintf(fd,"\n");
      }

      if (rs_nodes && (rs_node_cnt>0))
      {
         fprintf(fd,"nodalData = \n");
         fprintf(fd,"#  Primary lines: node_tag, x, y, z, num_arms, constraint\n");
         fprintf(fd,"#  Secondary lines: arm_tag, burgx, burgy, burgz, nx, ny, nz\n");

         Node_t *node = rs_nodes;

         for (int i=0; (i<rs_node_cnt); i++, node++)
         {
            fprintf(fd, " %4d,%-5d %14.8f %14.8f %14.8f %2d %2d\n",
                    node->myTag.domainID, node->myTag.index,
                    node->x, node->y, node->z, node->numNbrs, node->constraint);

            for (int j=0; (j<node->numNbrs); j++)
            {
               fprintf(fd,"    %4d,%-5d"                            , node->nbrTag[j].domainID, node->nbrTag[j].index );
               fprintf(fd," %11.8lf %11.8lf %11.8lf\n"              , node->burgX [j], node->burgY[j], node->burgZ[j] );
               fprintf(fd,"               %11.8lf %11.8lf %11.8lf\n", node->nx    [j], node->ny   [j], node->nz   [j] );
            }
         }
      }

   }
}

//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::Print (const char *path)
{
   FILE *fd = ( path ? fopen(path,"w") : 0 );

   if (fd) { Print(fd); fclose(fd); }
}

// Print_GNU
//
// Saves a gnuplot-compatible version of the current restart file.
//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::Print_GNU (const char *path)
{
   Print_GNU_Nodes(path);
   Print_GNU_Doms (path);
   Print_GNU_Cmd  (path);
}

void Paradis_Restart::Print_GNU_Doms (const char *path)
{
   char path_dom[256]; sprintf(path_dom, "%s.dom", path);

   if (rs_domains && (rs_domain_cnt>0))
   {
      FILE *fd = ( path ? fopen(path_dom,"w") : 0 );

      if (fd)
      {
         for (int i=0; (i<rs_domain_cnt); i++)
            rs_domains[i].Print_GNU(fd);

         fclose(fd);
      }
   }
}

void Paradis_Restart::Print_GNU_Cmd (const char *path)
{
   char path_cmd[256]; sprintf(path_cmd, "%s.gnu", path);

   FILE *fd = ( path ? fopen(path_cmd,"w") : 0 );

   int x0 = (int) rs_min_coords[0];  x0 = 1000*((x0/1000)-1);
   int y0 = (int) rs_min_coords[1];  y0 = 1000*((y0/1000)-1);
   int z0 = (int) rs_min_coords[2];  z0 = 1000*((z0/1000)-1);

   int x1 = (int) rs_max_coords[0];  x1 = 1000*((x1/1000)+1);
   int y1 = (int) rs_max_coords[1];  y1 = 1000*((y1/1000)+1);
   int z1 = (int) rs_max_coords[2];  z1 = 1000*((z1/1000)+1);

   if (fd)
   {
      fprintf(fd,"#set term x11\n"  );
      fprintf(fd,"#set term pngcairo size 1000,800\n");
      fprintf(fd,"#set out '%s.png'\n", path );
      fprintf(fd,"reset\n" );
      fprintf(fd,"set view 60,30\n");
      fprintf(fd,"set xrange [%d:%d]\n", x0, x1 );
      fprintf(fd,"set yrange [%d:%d]\n", y0, y1 );
      fprintf(fd,"set zrange [%d:%d]\n", z0, z1 );
      fprintf(fd,"unset key\n"     );
      fprintf(fd,"set title '%s.dat'\n",path);
      fprintf(fd,"set xlabel 'X'  \n" );
      fprintf(fd,"set ylabel 'Y'  \n" );
      fprintf(fd,"set zlabel 'Z'  \n" );
      fprintf(fd,"set ticslevel 0 \n" );
      fprintf(fd,"set xtics autofreq %d, %d, %d\n", x0, 1000, x1 );
      fprintf(fd,"set ytics autofreq %d, %d, %d\n", y0, 1000, y1 );
      fprintf(fd,"set ztics autofreq %d, %d, %d\n", z0, 1000, z1 );
      fprintf(fd,"z0=%d\n", z0 );
      fprintf(fd,"splot '%s.dat' u 5:6:7:11:12:13     w vec nohead lt 1 lw 1 lc rgb '#0000ff' , \\\n", path );
      fprintf(fd,"      '%s.dat' u 5:6:(z0):11:12:(0) w vec nohead lt 1 lw 1 lc rgb '#00aaaa' , \\\n", path );
      fprintf(fd,"      '%s.dat' u 5:6:7              w points pt 6 ps 0.5   lc rgb '#0000ff' , \\\n", path );
      fprintf(fd,"      '%s.dat' u 8:9:10             w points pt 6 ps 0.5   lc rgb '#0000ff' , \\\n", path );
      fprintf(fd,"      '%s.dom' u 2:3:4:5:6:7        w vec nohead lt 1 lw 1 lc rgb '#ff00ff' ; \n"  , path );
      fprintf(fd,"#exit\n"         );

      fclose(fd);
   }
}

void Paradis_Restart::Print_GNU_Nodes (const char *path)
{
   char path_dat[256]; sprintf(path_dat, "%s.dat", path);

   int        nseg = 0;
   Segment_t *segs = Segments(nseg);

   if (segs)
   {
      FILE *fd = ( path ? fopen(path_dat,"w") : 0 );

      if (fd)
      {
         for (int i=0; (i<nseg); i++)
         {
            Node_t *n1 = segs[i].node1;
            Node_t *n2 = segs[i].node2;
            real8  *p1 = segs[i].p1;
            real8  *p2 = segs[i].p2;

            real8   dx = (p2[0]-p1[0]);
            real8   dy = (p2[1]-p1[1]);
            real8   dz = (p2[2]-p1[2]);

            fprintf(fd,"%3d %-4d %3d %-4d %15.8lf %15.8lf %15.8lf %15.8lf %15.8lf %15.8lf  %15.8lf %15.8lf %15.8lf\n",
               n1->myTag.domainID, n1->myTag.index,
               n2->myTag.domainID, n2->myTag.index,
               p1[0],p1[1],p1[2],
               p2[0],p2[1],p2[2],
               dx, dy, dz );
         }

         fclose(fd);
      }
   }

   VDELETE(segs);
}

// Node_Count()
//
// Will return the number of nodes in a specified domain.
//---------------------------------------------------------------------------------------------------------

int Paradis_Restart::Node_Count (const int dndx)
{
   int ncnt=0;

   if (rs_nodes && (rs_node_cnt>0))
   {
      for (int i=0; (i<rs_node_cnt); i++)
         if (rs_nodes[i].myTag.domainID == dndx) { ncnt++; }
   }

   return(ncnt);
}

// Get_Nodes()
//
// Will return an array of nodes that occur within a given domain.
//---------------------------------------------------------------------------------------------------------

Node_t *Paradis_Restart::Get_Nodes (const int dndx, int & ncnt)
{
   Node_t *nodes = 0;
           ncnt  = 0;

   if (rs_nodes && (rs_node_cnt>0))
   {
      for (int i=0; (i<rs_node_cnt); i++)
         if (rs_nodes[i].myTag.domainID == dndx) { ncnt++; }

      nodes = ( (ncnt>0) ? new Node_t[ncnt] : 0 );

      if (nodes)
         for (int i=0,j=0; (i<rs_node_cnt); i++)
            if (rs_nodes[i].myTag.domainID == dndx) { nodes[j]=rs_nodes[i]; j++; }
   }

   return(nodes);
}

Node_t *Paradis_Restart::Get_Nodes (const Domain_t & dom, int & ncnt)
{
   Node_t *nodes = 0;
           ncnt  = 0;

   if (rs_nodes && (rs_node_cnt>0))
   {
      Node_t *p = rs_nodes;
      for (int i=0; (i<rs_node_cnt); i++, p++)
         if (dom.Inside(p->x,p->y,p->z) ) { ncnt++; }

      nodes = ( (ncnt>0) ? new Node_t[ncnt] : 0 );

      if (nodes)
         for (int i=0,j=0; (i<rs_node_cnt); i++)
            if (dom.Inside(p->x,p->y,p->z)) { nodes[j]=rs_nodes[i]; j++; }
   }

   return(nodes);
}

// Get_Nodes()
//
// Alternate version. Returns an array of node pointers where each entry is an array of nodes
// to nodes assigned to each domain.
//---------------------------------------------------------------------------------------------------------

Node_t **Paradis_Restart::Get_Nodes
(
   int **ncnts    ///< receives the resulting node counts, allocated and returned by this method
)
{
   Node_t **nodes=0;

   if (rs_nodes && (rs_node_cnt>0) && ncnts)
   {
      int n = Ndom();

      nodes = ( (n>0) ? new Node_t *[n] : 0 );
     *ncnts = ( (n>0) ? new int     [n] : 0 );

      if ( nodes) memset( nodes,0,n*sizeof(Node_t *));
      if (*ncnts) memset(*ncnts,0,n*sizeof(int     ));

      if (nodes && *ncnts)
      {
         for (int i=0; (i<n); i++)
            nodes[i] = Get_Nodes(i,(*ncnts)[i]);
      }
   }

   return(nodes);
}

// Index_Nodes()
//
// Will cycle through all the nodes in the restart file and give them a unique index (sets myIndex).
//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::Index_Nodes (void)
{
   if (rs_nodes && (rs_node_cnt>0) )
   {
      for (int i=0; (i<rs_node_cnt); i++) { rs_nodes[i].myIndex=i; }
   }
}

// Node_Array()
//
// Returns an array of node pointers to all the nodes within the restart file.
//---------------------------------------------------------------------------------------------------------

Node_t **Paradis_Restart::Node_Array (void)
{
   Node_t **narr = ( (rs_nodes && (rs_node_cnt>0)) ? new Node_t *[rs_node_cnt]: 0 );

   if (narr) { for (int i=0; (i<rs_node_cnt); i++) { narr[i] = rs_nodes+i; } }

   return(narr);
}

Node_t **Paradis_Restart::Node_Array (int & ncnt)
{
   Node_t **narr = ( (rs_nodes && (rs_node_cnt>0)) ? new Node_t *[rs_node_cnt]: 0 );
            ncnt = rs_node_cnt;

   if (narr) { for (int i=0; (i<rs_node_cnt); i++) { narr[i] = rs_nodes+i; } }

   return(narr);
}

// Node_Array()
//
// Returns an array of node pointers to the nodes within the restart file that have the specified
// domain index.
//---------------------------------------------------------------------------------------------------------

Node_t **Paradis_Restart::Node_Array (const int dndx, int & ncnt)
{
   Node_t **narr = 0;
            ncnt = 0;

   if (rs_nodes && (rs_node_cnt>0))
   {
      for (int i=0; (i<rs_node_cnt); i++)
         if (rs_nodes[i].myTag.domainID == dndx) { ncnt++; }

      narr = ( (ncnt>0) ? new Node_t *[ncnt] : 0 );

      if (narr)
         for (int i=0,j=0; (i<rs_node_cnt); i++)
            if (rs_nodes[i].myTag.domainID == dndx) { narr[j]=rs_nodes+i; j++; }
   }

   return(narr);
}

// ReIndex()
//
// Will re-index all the nodes/tags within the restart file. Note that all the new tags come from
// the specified domain.
//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::ReIndex (const int dndx)
{
   int    ntags    = rs_node_cnt;
   Tag_t *tags_old = ( (rs_nodes && (ntags>0)) ? new Tag_t[ntags] : 0 );
   Tag_t *tags_new = ( (rs_nodes && (ntags>0)) ? new Tag_t[ntags] : 0 );

   if (tags_old && tags_new)
   {
      for (int i=0; (i<ntags); i++) tags_old[i]=rs_nodes[i].myTag;
      for (int i=0; (i<ntags); i++) tags_new[i]=Tag_t(dndx,i);

      Node_t *node = rs_nodes;
      for (int i=0; (i<ntags); i++, node++)
      {
         node->myTag = tags_new[Tag_Index(tags_old,ntags,node->myTag)];

         for (int j=0; (j< node->numNbrs); j++)
            node->nbrTag[j] = tags_new[Tag_Index(tags_old,ntags,node->nbrTag[j])];
      }
   }

   VDELETE(tags_old);
   VDELETE(tags_new);
}

// Segment_Split()
//
// If an existing segment exceeds the maximum segment length, will
// split the segment by inserting a new node at the midpoint of the segment.
// When split - the original nodes are adjusted to point to the newly created
// split node.
//
// Returns NIL If no split occurs;
//
// Note - this routine is invoked by the Refine() method.  I've kept this function
// in the restart class because I'm making assumptions about how the newly created
// split node will be integrated into the existing node arrays. I'm also making
// assumptions on how the periodic boundaries will be handled.  This is not a generally
// useful topology function.  Care must be taken when appending newly created split nodes
// to the end of the node arrays because appending nodes can result in invalidating the
// node pointers associated with the segment list.
//
//   n1<------------>n2    (original   segment )
//   n1<---->n3<---->n2    (resulting  segments)
//---------------------------------------------------------------------------------------------------------

static Node_t *Segment_Split
(
   const Segment_t & seg ,   ///< the segment to be split
   const Tag_t     & tag ,   ///< tag for the new split node
   const real8       smax,   ///< maximum segment length
   const real8       lx  ,   ///< simulation size (x) (for PBC correction)
   const real8       ly  ,   ///< simulation size (y) (for PBC correction)
   const real8       lz      ///< simulation size (z) (for PBC correction)
)
{
   Node_t *node=0;

   Node_t *n1 = seg.node1;
   Node_t *n2 = seg.node2;

   if (n1 && n2)
   {
      real8 p1[3] = { n1->x, n1->y, n1->z };              // p1 = position of node 1
      real8 p2[3] = { n2->x, n2->y, n2->z };              // p2 = position of node 2

            p2[0] =  p2[0] - rint((p2[0]-p1[0])/lx)*lx;   // p2 = position of node 2 (corrected for PBC)
            p2[1] =  p2[1] - rint((p2[1]-p1[1])/ly)*ly;   //
            p2[2] =  p2[2] - rint((p2[2]-p1[2])/lz)*lz;   //

      real8 dx    = (p2[0]-p1[0]);
      real8 dy    = (p2[1]-p1[1]);
      real8 dz    = (p2[2]-p1[2]);

      // Only create a new segment if the segment length exceeds the max seg length...

      if ( sqrt(dx*dx + dy*dy + dz*dz) > smax )
      {
         p2[0] = (p1[0]+p2[0])/2.0;                  // p2 = midpoint of segment to split
         p2[1] = (p1[1]+p2[1])/2.0;                  //
         p2[2] = (p1[2]+p2[2])/2.0;                  //

         p2[0] = p2[0] - rint(p2[0]/lx)*lx;          // p2 = fold back into simulation box
         p2[1] = p2[1] - rint(p2[1]/ly)*ly;          //
         p2[2] = p2[2] - rint(p2[2]/lz)*lz;          //

         node  = new Node_t(tag,2);                  // create a new node with two arms using the provided tag

         if (node)
         {
            const Tag_t t1  = n1->myTag;             // t1 = tag of node 1
            const Tag_t t2  = n2->myTag;             // t2 = tag of node 2

            const int   j1  = n1->Arm_Index(t2);     // j1 = index of arm t2 of node 1
            const int   j2  = n2->Arm_Index(t1);     // j2 = index of arm t1 of node 2

            node->x         = p2[0];                 // set split node position
            node->y         = p2[1];                 //
            node->z         = p2[2];                 //

            node->nbrTag[0] = t1;                    // set split node arm tags
            node->nbrTag[1] = t2;                    //

            node->burgX [0] = n2->burgX[j2];         // set split node arm burgers vectors
            node->burgY [0] = n2->burgY[j2];         //
            node->burgZ [0] = n2->burgZ[j2];         //

            node->burgX [1] = n1->burgX[j1];         //
            node->burgY [1] = n1->burgY[j1];         //
            node->burgZ [1] = n1->burgZ[j1];         //

            node->nx    [0] = n2->nx   [j2];         // set split node arm normals
            node->ny    [0] = n2->ny   [j2];         //
            node->nz    [0] = n2->nz   [j2];         //

            node->nx    [1] = n1->nx   [j1];         //
            node->ny    [1] = n1->ny   [j1];         //
            node->nz    [1] = n1->nz   [j1];         //

            n1->nbrTag [j1] = tag;                   // update the existing nodes arms to the split node
            n2->nbrTag [j2] = tag;                   //
         }
      }
   }

   return(node);
}

// Refine()
//
// Given an existing configuration, will cycle through all the existing segments
// and bisect those that have segment lengths that exceed the provided maximum.
//
// Note - single domain only.
//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::Refine
(
   const real8 smax   ///< maximum segment length (b)
)
{
   ReIndex(0); // (start with a new, clean node index, domain zero)

   const real8 lx   = Lx();
   const real8 ly   = Ly();
   const real8 lz   = Lz();

         int   ncnt = rs_node_cnt;
         int   nmax = rs_node_cnt;

   // The strategy with each iteration is to build a segment list
   // from the existing node configuration. We then loop through
   // all the segments and split those segments that exceed the
   // maximum segment length provided. This will continue until there are
   // no segments remaining to split.

   for (int nsplit=1; (nsplit>0); )
   {
                 nsplit = 0;                 // reset split count for this iteration

      int        nseg   = 0;                 //
      Segment_t *segs   = Segments(nseg);    // build new segment list

      if (segs && (nseg>0))
      {
         // Note that we cannot simply append new nodes to the end of the current
         // node list because doing so could potentially invalidate all the node
         // pointers in the current segment list.  Instead, we will maintain a
         // separate list of split nodes and then append them after we've visited
         // all the segments.

         Node_t    *split_nodes = 0;         // local array of split nodes
         int        split_ncnt  = 0;         //
         int        split_nmax  = 0;         //

         for (int i=0; (i<nseg); i++)
         {
            Tag_t   tag  = Tag_t(0,ncnt);
            Node_t *node = Segment_Split(segs[i],tag,smax,lx,ly,lz);

            if (node)
            {
               split_nodes = Node_Append(split_nodes, split_ncnt, split_nmax, *node);
               ncnt++;       // increment total node count (used to bump the next tag index)
               nsplit++;     // increment split node count
            }
         }

         // append newly created split nodes to the main node list...

         if (split_nodes)
         {
            for (int i=0; (i<split_ncnt); i++)
            { rs_nodes = Node_Append(rs_nodes, rs_node_cnt, nmax, split_nodes[i]); }

            delete [] split_nodes; split_nodes=0; split_ncnt=0; split_nmax=0;
         }
      }

      if (segs) { delete [] segs; segs=0; nseg=0; }  // (cleanup)
   }

   Sort_Nodes  ();
   Index_Nodes ();
}

// Segments()
//
// Will construct and return an array of unique segments in the restart file.
// The method assumes the nodes have been sorted. A segment is considered unique when
// the tag of the fist node is less than the tag of the second node.
//---------------------------------------------------------------------------------------------------------

Segment_t  *Paradis_Restart::Segments (int & nseg)
{
   Segment_t *segs = 0;
              nseg = 0;

   const real8 lx = rs_max_coords[0]-rs_min_coords[0];
   const real8 ly = rs_max_coords[1]-rs_min_coords[1];
   const real8 lz = rs_max_coords[2]-rs_min_coords[2];

   if (rs_nodes && (rs_node_cnt>0))
   {
      // first pass...determine the number of segments
      // note - during this pass, we are only counting the number
      // of possible segments.

      Node_t *n1 = (Node_t *) rs_nodes;

      for (int i=0; (i<rs_node_cnt); i++, n1++)
      {
         for (int j=0; (j< n1->numNbrs); j++)
            if (n1->myTag < n1->nbrTag[j]) { nseg++; }
      }

      // allocate space for the segments...

      segs = ( (nseg>0) ? new Segment_t[nseg] : 0 );

      // second pass...build the atual segment array...
      // note - during this pass we actually perform the
      // search for the second node. It's possible that more
      // segments were allocated than actually found.

      nseg = 0; // (reset to get actual count during this pass)

      n1 = (Node_t *) rs_nodes;

      for (int i=0; (i<rs_node_cnt); i++, n1++)
      {
         for (int j=0; (j< n1->numNbrs); j++)
         {
            if (n1->myTag < n1->nbrTag[j])
            {
               Node_t *n2 = Find_Node(n1->nbrTag[j]);

               if (n2)
               {
                  real8 bv[3] = { n1->burgX[j], n1->burgY[j], n1->burgZ[j] };  // bv = burgers vector
                  real8 nv[3] = { n1->nx   [j], n1->ny   [j], n1->nz   [j] };  // nv = glide plane

                  segs[nseg] = Segment_t(n1,n2,bv,nv);
                  segs[nseg].PBC_Position(lx,ly,lz);  nseg++;
               }
            }
         }
      }
   }

   return(segs);
}

// Segments()
//
// Will extract a subset of the unique segments specific to a domain.
//---------------------------------------------------------------------------------------------------------

Segment_t  *Paradis_Restart::Segments (int & nseg, const int dndx)
{
   Segment_t *segs = 0;
              nseg = 0;

   const real8 lx = rs_max_coords[0]-rs_min_coords[0];
   const real8 ly = rs_max_coords[1]-rs_min_coords[1];
   const real8 lz = rs_max_coords[2]-rs_min_coords[2];

   if (rs_nodes && (rs_node_cnt>0))
   {
      // first pass...determine the number of segments
      // note - during this pass, we are only counting the number
      // of possible segments.

      Node_t *n1 = (Node_t *) rs_nodes;

      for (int i=0; (i<rs_node_cnt); i++, n1++)
      {
         for (int j=0; (j< n1->numNbrs); j++)
            if (    (n1->myTag < n1->nbrTag[j])
                 && ( (n1->myTag.Domain() == dndx) || (n1->nbrTag[j].Domain() == dndx) ) )
            { nseg++; }
      }

      // allocate space for the segments...

      segs = ( (nseg>0) ? new Segment_t[nseg] : 0 );

      // second pass...build the actual segment array...
      // note - during this pass we actually perform the
      // search for the second node. It's possible that more
      // segments were allocated than actually found.

      nseg = 0; // (reset to get actual count during this pass)

      n1 = (Node_t *) rs_nodes;

      for (int i=0; (i<rs_node_cnt); i++, n1++)
      {
         for (int j=0; (j< n1->numNbrs); j++)
         {
            if (    (n1->myTag < n1->nbrTag[j])
                 && ( (n1->myTag.Domain() == dndx) || (n1->nbrTag[j].Domain() == dndx) ) )
            {
               Node_t *n2 = Find_Node(n1->nbrTag[j]);

               if (n2)
               {
                  real8 bv[3] = { n1->burgX[j], n1->burgY[j], n1->burgZ[j] };  // bv = burgers vector
                  real8 nv[3] = { n1->nx   [j], n1->ny   [j], n1->nz   [j] };  // nv = glide plane

                  segs[nseg] = Segment_t(n1,n2,bv,nv);
                  segs[nseg].PBC_Position(lx,ly,lz);  nseg++;
               }
            }
         }
      }
   }

   return(segs);
}

// Validate()
//
// This is a quick validation diagnostic.
//
// Returns (1==valid) if the sum of the burgers vectors for all the node arms is zero.
//---------------------------------------------------------------------------------------------------------

int Paradis_Restart::Validate (const real8 eps)
{
   if (rs_nodes && (rs_node_cnt>0))
   {
      Node_t *n = rs_nodes;

      for (int i=0; (i<rs_node_cnt); i++, n++)
         if (!n->Validate()) return(0);
   }

   return(1);
}

//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::Reset_Forces (void)
{
   if (rs_nodes && (rs_node_cnt>0))
   {
      for (int i=0; (i<rs_node_cnt); i++) { rs_nodes[i].Forces_Reset(); }
   }
}

// Max_Arms()
//
// Returns the maximum arm count of all the nodes (used by the arm histogram method).
//---------------------------------------------------------------------------------------------------------

int  Paradis_Restart::Max_Arms (void)
{
   int max=0;

   if (rs_nodes && (rs_node_cnt>0))
   {
      for (int i=0; (i<rs_node_cnt); i++)
         max = ( (rs_nodes[i].numNbrs>max) ? rs_nodes[i].numNbrs : max );
   }

   return(max);
}

// Arm_Histogram()
//
// Returns a histogram of all the arms in the dataset.
//---------------------------------------------------------------------------------------------------------

int *Paradis_Restart::Arm_Histogram (int & max_arms)
{
   max_arms = Max_Arms();

   int *histo = ( (max_arms>0) ? new int[max_arms+1] : 0 );

   if (histo)
   {
      memset(histo,0,(max_arms+1)*sizeof(int));

      for (int i=0; (i<rs_node_cnt); i++)
      {
         int narms = rs_nodes[i].numNbrs;

         if ( (0<=narms) && (narms<=max_arms) ) { histo[narms]++; }
      }
   }

   return(histo);
}

//---------------------------------------------------------------------------------------------------------
// Burgers vector table generation.
//---------------------------------------------------------------------------------------------------------

// BV_Index()
//
// Given an array of burgers vectors, will return the index into the array of the desired burgers
// vector.  Returns -1 if the vector is not in the list.
//---------------------------------------------------------------------------------------------------------

static int BV_Index
(
         real8  *bv,        ///< current burgers vector table
         int     bn,        ///< current burgers vector count
   const real8   bx,        ///< burgers vector <x>
   const real8   by,        ///< burgers vector <y>
   const real8   bz,        ///< burgers vector <z>
   const real8   eps=1.0e-6 ///< error tolerance for comparison
)
{
   if (bv && (bn>0))
   {
      for (int i=0,j=0; (i<bn); i++, j+=3)
      {
         if (    (fabs(bv[j  ]-bx)<eps)
              && (fabs(bv[j+1]-by)<eps)
              && (fabs(bv[j+2]-bz)<eps) ) return(i);
      }
   }

   return(-1);  // (not found)
}

// BV_Append()
//
// Will append a new burgers vector to an existing array.  Manages the allocation as necessary.
//---------------------------------------------------------------------------------------------------------

static real8 *BV_Append
(
         real8  *bv  ,   ///< current burgers vector table
         int   & bn  ,   ///< current burgers vector count
         int   & bmax,   ///< length of current allocation
   const real8   bx  ,   ///< burgers vector <x>
   const real8   by  ,   ///< burgers vector <y>
   const real8   bz      ///< burgers vector <z>
)
{
   if ( !bv || (bn==bmax) )
   {
      bmax = ( (bmax==0) ? 32 : 2*bmax );

      real8 *tmp = new real8[3*bmax];

      if (tmp && bv) { memcpy(tmp,bv,3*bn*sizeof(real8)); }

      VDELETE(bv);

      bv = tmp;
   }

   if (bv && (bn<bmax))
   {
      int j=3*bn; bn++;
      bv[j  ] = bx;
      bv[j+1] = by;
      bv[j+2] = bz;
   }

   return(bv);
}

// Burgers_Vectors()
//
// Will return a pointer to a list of all the burgers vectors encountered in the restart file.
//---------------------------------------------------------------------------------------------------------

real8 *Paradis_Restart::Burgers_Vectors (void)
{
   real8 *bv=0; int bmax=0, bcnt=0;

   if (rs_nodes && (rs_node_cnt>0))
   {
      for (int i=0; (i<rs_node_cnt); i++)
      {
         int narms = rs_nodes[i].numNbrs;

         for (int j=0; (j<narms); j++)
         {
            real8 bx = rs_nodes[i].burgX[j];
            real8 by = rs_nodes[i].burgY[j];
            real8 bz = rs_nodes[i].burgZ[j];

            // check and add both bv and -bv to the table...

            if ( BV_Index(bv,bcnt, bx, by, bz)<0)  { bv = BV_Append(bv,bcnt,bmax, bx, by, bz); }
            if ( BV_Index(bv,bcnt,-bx,-by,-bz)<0)  { bv = BV_Append(bv,bcnt,bmax,-bx,-by,-bz); }
         }
      }
   }

   VDELETE(rs_bv_tbl);

   rs_bv_cnt = bcnt;
   rs_bv_tbl = bv;

   return(bv);
}

// Burgers_Lengths()
//
// Returns the total length of all segments within the restart file, sorted by burgers vector.
//---------------------------------------------------------------------------------------------------------

real8 *Paradis_Restart::Burgers_Lengths (void)
{
   real8      *bv    = ( rs_bv_tbl ? rs_bv_tbl : Burgers_Vectors() );
   int         bcnt  = rs_bv_cnt;
   real8      *blens = ( (bcnt>0) ? new real8[bcnt] : 0 );

   int         scnt  = 0;
   Segment_t  *segs  = Segments(scnt);

   const real8 lx = Lx();
   const real8 ly = Ly();
   const real8 lz = Lz();

   const real8 sx = ( (fabs(lx)>0.0) ? 1.0/fabs(lx) : 0.0 );
   const real8 sy = ( (fabs(ly)>0.0) ? 1.0/fabs(ly) : 0.0 );
   const real8 sz = ( (fabs(lz)>0.0) ? 1.0/fabs(lz) : 0.0 );

   if (blens) memset(blens,0,bcnt*sizeof(real8));

   if (bv && blens && segs && (scnt>0))
   {
      for (int i=0; (i<scnt); ++i)
      {
         real8 bx = segs[i].bv[0];
         real8 by = segs[i].bv[1];
         real8 bz = segs[i].bv[2];

         int         bndx = BV_Index(bv,bcnt, bx, by, bz);
         if (bndx<0) bndx = BV_Index(bv,bcnt,-bx,-by,-bz);

         if (bndx>=0)
         {
            real8 p1[3] = { segs[i].p1[0], segs[i].p1[1], segs[i].p1[2] };
            real8 p2[3] = { segs[i].p2[0], segs[i].p2[1], segs[i].p2[2] };

            PBC_Position(p1,p2, lx,ly,lz, sx,sy,sz);

            real8 dx = (p2[0]-p1[0]);
            real8 dy = (p2[1]-p1[1]);
            real8 dz = (p2[2]-p1[2]);

            blens[bndx] += sqrt( (dx*dx) + (dy*dy) + (dz*dz) );
         }
      }
   }

   VDELETE(segs);

   return (blens);
}

//---------------------------------------------------------------------------------------------------------

void Paradis_Restart::Extract
(
   const Paradis_Restart & rs,   ///< source restart file
   const real8 *p ,              ///< center of the extraction window
   const real8 *pw               ///< size of the extraction window
)
{
   // start with a clean slate...

   Recycle();

   rs_version = 5;
   rs_nsegs   = 1;

   // create the physical window...

   const real8 x0 = (p[0]-pw[0]/2.0);
   const real8 x1 = (p[0]+pw[0]/2.0);
   const real8 y0 = (p[1]-pw[1]/2.0);
   const real8 y1 = (p[1]+pw[1]/2.0);
   const real8 z0 = (p[2]-pw[2]/2.0);
   const real8 z1 = (p[2]+pw[2]/2.0);

   rs_min_coords[0]=x0;
   rs_min_coords[1]=y0;
   rs_min_coords[2]=z0;

   rs_max_coords[0]=x1;
   rs_max_coords[1]=y1;
   rs_max_coords[2]=z1;

   rs_domains    = new Domain_t[1];
   rs_domain_cnt = 1;

   rs_domains[0] = Domain_t(0, x0,x1, y0,y1, z0,z1);

   rs_decomp_type   = 0;
   rs_decomp_geo[0] = 1;
   rs_decomp_geo[1] = 1;
   rs_decomp_geo[2] = 1;

   // identify all nodes within the physical window...

   int nmax=0;

   for (int i=0; (i<rs.rs_node_cnt); i++)
   {
      if ( rs.rs_nodes[i].Inside(x0,x1, y0,y1, z0,z1) )
      { rs_nodes = Node_Append(rs_nodes, rs_node_cnt, nmax, rs.rs_nodes[i]); }
   }

   // identify all the ghost nodes that reside outside the window...

   Node_t **narr0 = Nodes_Vector(rs.rs_nodes, rs.rs_node_cnt);    // narr0 = array of node pointers to full node list
   Node_t **narr1 = Nodes_Vector(   rs_nodes,    rs_node_cnt);    // narr1 = array of node pointers to extracted list

   Nodes_Sort(narr0,rs.rs_node_cnt);  // (sort nodes for lookups)
   Nodes_Sort(narr1,   rs_node_cnt);  //

   // build up a an array of nodes that are outside of the extraction window, but connected to nodes
   // within the window

   Node_t *ghosts=0; int gcnt=0, gmax=0;

   if (narr0 && narr1 && rs_nodes && rs.rs_nodes && (rs_node_cnt>0) )
   {
      // loop through all the extracted nodes...

      for (int i=0; (i<rs_node_cnt); i++)
      {
         int    narms = rs_nodes[i].numNbrs;
         Tag_t *tags  = rs_nodes[i].nbrTag;

         if ( tags && (narms>0) )
         {
            // visit all the arms of the extracted nodes...

            for (int j=0; (j<narms); j++)
            {
               // if the arm node isn't in the current list, append to the ghosts...

               if ( !Nodes_Find(narr1,rs_node_cnt,tags[j]) )
               {
                  Node_t *pn = Nodes_Find(narr0,rs.rs_node_cnt,tags[j]);

                  if (pn) { ghosts=Node_Append(ghosts,gcnt,gmax,*pn); }
               }
            }
         }
      }
   }

   // append all the ghosts to the extracted node_array...

   if (ghosts && (gcnt>0))
   {
      for (int i=0; (i<gcnt); i++)
      { rs_nodes=Node_Append(rs_nodes,rs_node_cnt,nmax,ghosts[i]); }
   }

   // We now have a complete node list that includes the nodes within the extraction region
   // and any nodes outside the region that are connected to those nodes.

   Sort_Nodes();  // (sort the updated list for lookups)

   // Scan through all the nodes and delete any arms that are not connected to nodes
   // within the extraction window.

   if (rs_nodes && (rs_node_cnt>0))
   {
      for (int i=0; (i<rs_node_cnt); i++)
      {
         int    narms = rs_nodes[i].numNbrs;
         Tag_t *tags  = rs_nodes[i].nbrTag;

         if ( tags && (narms>0) )
         {
            for (int j=0; (j<narms); j++)
            {
               if ( !Nodes_Find(rs_nodes,rs_node_cnt,tags[j]) )
               {
                  rs_nodes[i].Delete_Arm(tags[j]);
                  rs_nodes[i].constraint = PINNED_NODE;
                  narms--; j--;
               }
            }
         }
      }
   }

   // cleanup...

   VDELETE(narr0 );
   VDELETE(narr1 );
   VDELETE(ghosts);  gcnt=0; gmax=0;
}
