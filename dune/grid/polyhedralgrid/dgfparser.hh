// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_DGFPARSER_HH
#define DUNE_POLYHEDRALGRID_DGFPARSER_HH

#include <algorithm>
#include <numeric>

#include <dune/common/typetraits.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/grid/polyhedralgrid/grid.hh>

namespace Dune
{

  namespace dgf
  {

    namespace PolyhedralGrid
    {

      // PolygonBlock
      // ------------

      struct PolygonBlock
        : public BasicBlock
      {
        PolygonBlock ( std::istream &in, int numVtx, int vtxOfs )
          : BasicBlock( in, "Polygon" ), vtxBegin_( vtxOfs ), vtxEnd_( vtxOfs + numVtx )
        {}

        int get ( std::vector< std::vector< int > > &polygons )
        {
          reset();
          std::vector< int > polygon;
          while( getnextline() )
          {
            polygon.clear();
            for( int vtxIdx; getnextentry( vtxIdx ); )
            {
              if( (vtxBegin_ > vtxIdx) || (vtxIdx >= vtxEnd_) )
                DUNE_THROW( DGFException, "Error in " << *this << ": Invalid vertex index (" << vtxIdx << " not int [" << vtxBegin_ << ", " << vtxEnd_ << "[)" );
              polygon.push_back( vtxIdx - vtxBegin_ );
            }

            polygons.push_back( polygon );
          }
          return polygons.size();
        }

      private:
        int vtxBegin_, vtxEnd_;
      };



      // PolyhedronBlock
      // ---------------

      struct PolyhedronBlock
        : public BasicBlock
      {
        explicit PolyhedronBlock ( std::istream &in, int numPolys )
          : BasicBlock( in, "Polyhedron" ), numPolys_( numPolys )
        {}

        int get ( std::vector< std::vector< int > > &polyhedra )
        {
          reset();
          std::vector< int > polyhedron;
          int minPolyId = 1;
          while( getnextline() )
          {
            polyhedron.clear();
            for( int polyIdx; getnextentry( polyIdx ); )
            {
              if( (polyIdx < 0) || (polyIdx > numPolys_) )
                DUNE_THROW( DGFException, "Error in " << *this << ": Invalid polygon index (" << polyIdx << " not int [0, " << numPolys_ << "])" );

              minPolyId = std::min( minPolyId, polyIdx );
              polyhedron.push_back( polyIdx );
            }

            polyhedra.push_back( polyhedron );
          }

          // substract minimal number to have 0 starting numbering
          if( minPolyId > 0 )
          {
            const size_t polySize = polyhedra.size();
            for( size_t i=0; i<polySize; ++i )
            {
              const size_t pSize = polyhedra[ i ].size();
              for( size_t j=0; j<pSize; ++j )
              {
                polyhedra[ i ][ j ] -= minPolyId;
              }
            }
          }
          return polyhedra.size();
        }

      private:
        const int numPolys_;
      };

    } // namespace PolyhedralGrid

  } // namespace dgf



  // DGFGridFactory for PolyhedralGrid
  // ---------------------------------

  template< int dim, int dimworld >
  struct DGFGridFactory< PolyhedralGrid< dim, dimworld > >
  {
    typedef PolyhedralGrid< dim, dimworld > Grid;

    const static int dimension = Grid::dimension;
    typedef MPIHelper::MPICommunicator MPICommunicator;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;

    explicit DGFGridFactory ( std::istream &input, MPICommunicator comm = MPIHelper::getCommunicator() )
      : grid_( nullptr )
    {
      input.clear();
      input.seekg( 0 );
      if( !input )
        DUNE_THROW( DGFException, "Error resetting input stream" );
      generate( input );
    }

    explicit DGFGridFactory ( const std::string &filename, MPICommunicator comm = MPIHelper::getCommunicator() )
      : grid_( nullptr )
    {
      std::ifstream input( filename );
      if( !input )
        DUNE_THROW( DGFException, "Macrofile '" << filename << "' not found" );
      generate( input );
    }

    Grid *grid () const { return grid_; }

    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return false;
    }

    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      return false;
    }

    bool haveBoundaryParameters () const { return false; }

    template< int codim >
    int numParameters () const
    {
      //return (codim == dimension ? numVtxParams_ : 0);;
      return 0;
    }

    template< class Intersection >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection &intersection ) const
    {
      return DGFBoundaryParameter::defaultValue();;
    }

    template< class Entity >
    std::vector< double > &parameter ( const Entity &entity )
    {
      static std::vector< double > dummy;
      return dummy;
    }

  private:
    int readVertices ( std::istream &input, std::vector< std::vector< double > > &vertices )
    {
      int dimWorld = 3;
      dgf::VertexBlock vtxBlock( input, dimWorld );
      if( !vtxBlock.isactive() )
        DUNE_THROW( DGFException, "Vertex block not found" );

      vtxBlock.get( vertices, vtxParams_, numVtxParams_ );
      return vtxBlock.offset();
    }

    std::vector< std::vector< int > > readPolygons ( std::istream &input, int numVtx, int vtxOfs )
    {
      dgf::PolyhedralGrid::PolygonBlock polygonBlock( input, numVtx, vtxOfs );
      if( !polygonBlock.isactive() )
        DUNE_THROW( DGFException, "Polygon block not found" );

      std::vector< std::vector< int > > polygons;
      polygonBlock.get( polygons );
      return polygons;
    }

    std::vector< std::vector< int > > readPolyhedra ( std::istream &input, int numPolygons )
    {
      dgf::PolyhedralGrid::PolyhedronBlock polyhedronBlock( input, numPolygons );
      if( !polyhedronBlock.isactive() )
        DUNE_THROW( DGFException, "Polyhedron block not found" );

      std::vector< std::vector< int > > polyhedra;
      polyhedronBlock.get( polyhedra );
      return polyhedra;
    }

    template< class Iterator >
    void copy ( Iterator begin, Iterator end, double *dest )
    {
      for( ; begin != end; ++begin )
        dest = std::copy( begin->begin(), begin->end(), dest );
    }

    template< class Iterator >
    void copy ( Iterator begin, Iterator end, int *dest, int *offset )
    {
      int size = 0;
      for( ; begin != end; ++begin )
      {
        *(offset++) = size;
        size += begin->size();
        dest = std::copy( begin->begin(), begin->end(), dest );
      }
      *offset = size;
    }

    void generate ( std::istream &input )
    {
      if( !DuneGridFormatParser::isDuneGridFormat( input ) )
        DUNE_THROW( DGFException, "Not in DGF format" );

      std::vector< std::vector< double > > nodes;
      const int vtxOfs = readVertices( input, nodes );

      std::vector< std::vector< int > > faces = readPolygons( input, nodes.size(), vtxOfs );
      std::vector< std::vector< int > > cells = readPolyhedra( input, faces.size() );

      const auto sumSize = [] ( std::size_t s, const std::vector< int > &v ) { return s + v.size(); };
      const std::size_t numFaceNodes = std::accumulate( faces.begin(), faces.end(), std::size_t( 0 ), sumSize );
      const std::size_t numCellFaces = std::accumulate( cells.begin(), cells.end(), std::size_t( 0 ), sumSize );

      typename Grid::UnstructuredGridPtr ug = Grid::allocateGrid( cells.size(), faces.size(), numFaceNodes, numCellFaces, nodes.size() );

      // copy faces
      {
#ifndef NDEBUG
        std::map< std::vector< int >, std::vector< int > > faceMap;
#endif

        const int nFaces = faces.size();
        // set all face_cells values to -2 as default
        std::fill( ug->face_cells, ug->face_cells + 2*nFaces, -1 );

        int facepos = 0;
        std::vector< int > faceVertices;
        faceVertices.reserve( 30 );
        for( int face = 0; face < nFaces; ++face )
        {
          //std::cout << "face " << face << ": ";
          faceVertices.clear();
          ug->face_nodepos[ face ] = facepos;
          const int nVertices = faces[ face ].size();
          for( int vx = 0; vx < nVertices; ++vx, ++facepos )
          {
            //std::cout << " " << faces[ face ][ vx ];
            ug->face_nodes[ facepos ] = faces[ face ][ vx ];
            faceVertices.push_back( faces[ face ][ vx ] );
          }
          //std::cout << std::endl;

#ifndef NDEBUG
          // sort vertices
          std::sort( faceVertices.begin(), faceVertices.end() );
          // make sure each face only exists once
          faceMap[ faceVertices ].push_back( face );
          assert( faceMap[ faceVertices ].size() == 1 );
#endif
        }
        ug->face_nodepos[ nFaces ] = facepos ;
      }

      // copy cells
      {
        const int nCells = cells.size();
        int cellpos = 0;
        for( int cell = 0; cell < nCells; ++cell )
        {
          //std::cout << "Cell " << cell << ": ";
          ug->cell_facepos[ cell ] = cellpos;
          const int nFaces = cells[ cell ].size();
          for( int f = 0; f < nFaces; ++f, ++cellpos )
          {
            const int face = cells[ cell ][ f ];
            // std::cout << " " << face ;
            ug->cell_faces[ cellpos ] = face;

            // TODO find cells for each face
            if( ug->face_cells[ 2*face ] == -1 )
            {
              ug->face_cells[ 2*face ] = cell;
            }
            else // if ( ug->face_cells[ 2*face+1 ] == -1 )
            {
              //assert( ug->face_cells[ 2*face+1 ] == -1 );
              ug->face_cells[ 2*face+1 ] = cell;
            }
          }
          //std::cout << std::endl;
        }
        ug->cell_facepos[ nCells ] = cellpos ;
      }

      // copy node coordinates
      {
        const int nNodes = nodes.size();
        int nodepos = 0;
        for( int vx = 0 ; vx < nNodes; ++vx )
        {
          for( int d=0; d<dim; ++d, ++nodepos )
            ug->node_coordinates[ nodepos ] = nodes[ vx ][ d ];
        }
      }

      /*
      for( int i=0; i<int(faces.size() ); ++i)
      {
        std::cout << "face "<< i<< " connects to " << ug->face_cells[ 2*i ] << " " <<
          ug->face_cells[ 2*i+1] << std::endl;
      }
      */

      // free cell face tag since it's not a cartesian grid
      if( ug->cell_facetag )
      {
        std::free( ug->cell_facetag );
        ug->cell_facetag = nullptr ;
        for( int i=0; i<3; ++i ) ug->cartdims[ i ] = 0;
      }

      // compute geometric quntities like cell volume and face normals
      Grid::computeGeometry( ug );

      grid_ = new Grid( std::move( ug ) );
    }

    Grid *grid_;
    int numVtxParams_;
    std::vector< std::vector< double > > vtxParams_;
  };



  // DGFGridInfo for PolyhedralGrid
  // ------------------------------

  template< int dim, int dimworld >
  struct DGFGridInfo< PolyhedralGrid< dim, dimworld > >
  {
    static int refineStepsForHalf ()
    {
      return 0;
    }

    static double refineWeight ()
    {
      return 0;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_DGFPARSER_HH
