/*
  Copyright 2014 IRIS AS 

  This file is part of The Open Porous Media project  (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/CpGrid.hpp>
#include <opm/core/io/eclipse/EclipseGridInspector.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/common/fromtogridfactory.hh>

typedef Dune::ALUGrid< 3, 3, Dune::cube, Dune::nonconforming > GridType;

GridType* deck2dune( Opm::DeckConstPtr deck ) 
{
    // create empty CpGrid
    Dune::CpGrid cpgrid;

    // create CpGrid from deck
    cpgrid.processEclipseFormat(deck, 0.0, false, false);

    // grid factory converting a grid
    Dune::FromToGridFactory< GridType > factory;

    // create Grid from CpGrid
    Dune::GridPtr< GridType > grid = factory.convert( cpgrid );

    // return pointer
    return grid.release();
}

int main(int argc, char** argv)
try
{

    Dune::MPIHelper::instance( argc, argv );

    if (argc != 2) {
        std::cout << "Usage: grdecl2vtu filename.grdecl" << std::endl;
        exit(1);
    }
    
    const char* eclipsefilename = argv[1];
    Opm::ParserPtr parser(new Opm::Parser());
    Opm::DeckConstPtr deck(parser->parseFile(eclipsefilename));

    // create DUNE grid
    Dune::GridPtr< GridType > gridptr = Dune::GridPtr< GridType >(deck2dune( deck ));
    GridType& grid = *gridptr;
    
    Dune::VTKWriter<GridType::LeafGridView> vtkwriter(grid.leafGridView());

    std::string fname(eclipsefilename);
    std::string fnamebase = fname.substr(0, fname.find_last_of('.'));
    std::cout << "Writing to filename " << fnamebase << ".vtu" << std::endl;
    vtkwriter.write(fnamebase, Dune::VTK::ascii);
}
catch (const std::exception &e) {
    std::cerr << "Program threw an exception: " << e.what() << "\n";
    throw;
}
#else
int main() { return 0; }
#endif

