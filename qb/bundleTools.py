# Copyright (C) 2019  Andrea Vázquez Varela

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


#from soma import aims
import numpy as N
import os
import shutil
import math
#from soma import aims


def read_bundle( infile ):

  points = []
  bun_file = infile + 'data'
  os.path.getsize( bun_file )
  bytes = os.path.getsize( bun_file )
  num = bytes / 4

  num_count = 0
 # f = open( bun_file )
  with open(bun_file,'rb') as f:
    while num_count < num:
      # numero de puntos de la fibra

      p = N.frombuffer(f.read(4), 'i' )[ 0 ]
      vertex = N.frombuffer( f.read( p * 3 * 4 ), 'f' ).reshape( -1, 3 ) # lee coordenadas fibra
      points.append( vertex )
      #num_# commentunt = num_count + 1 + ( p * 3 )
      num_count = num_count + 1 + ( p * 3 )
#  f.close()

  #bundles = []
  return points#, bundles


def write_bundle( outfile, points ):

  #write bundles file
  f = open( outfile + 'data','wb' )
  ncount = len( points )
  for i in range( ncount ):
    f.write(N.array( [ len( points[ i ] ) ], N.int32 ).tostring() )
    f.write( points[ i ].ravel().tostring() )

  f.close()

  # wrtie minf file
  minf = """attributes = {\n    'binary' : 1,\n    'bundles' : %s,\n    'byte_order' : 'DCBA',\n    'curves_count' : %s,\n    'data_file_name' : '*.bundlesdata',\n    'format' : 'bundles_1.0',\n    'space_dimension' : 3\n  }"""
  open( outfile, 'w' ).write(minf % ( [ 'points', 0 ], ncount ) )



    
   
