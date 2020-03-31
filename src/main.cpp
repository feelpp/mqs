/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil;
 c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2014-06-03

 Copyright (C) 2014-2019 Cemosis - University of Strasbourg

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include <feel/feelcore/environment.hpp>

int main( int argc, char **argv )
{
  //# marker1 #
  using namespace Feel;
  po::options_description mqsoptions( "MQS options" );
  mqsoptions.add_options()
      ( "case.dim", po::value<int>()->default_value( 2 ), "dimension of the case study" );

  Environment env( _argc = argc, _argv = argv, _desc = mqsoptions,
                   _about = about( _name = "mqs",
                                   _author = "Feel++ Consortium",
                                   _email = "feelpp-devel@feelpp.org" ) );

}