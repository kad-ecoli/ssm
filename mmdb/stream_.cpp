//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2008.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    29.01.10   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  CStream_ <interface>
//       ~~~~~~~~~
//  **** Classes :  CStream  ( Basic Stream Class )
//       ~~~~~~~~~
//
//   (C) E. Krissinel 1995-2010
//
//  =================================================================
//

#ifndef  __Stream__
#include "stream_.h"
#endif

//  ==========================  CStream  ===========================

//     Each streamable class should be derived from CStream
//  and have constructor CClass(PCStream & Object), which should
//  initialize all memory of the class, and virtual functions
//  read(..) and write(..) (see below). Constructor CClass(PCStream&)
//  must not touch the Object variable. This constructor is used
//  only once just before read(..) function. It is assumed that
//  read/write functions of CClass provide storage/reading of
//  all vital data. Function read(..) must read data in exactly
//  the same way as function write(..) stores it.
//     For using CClass in streams, three following functions should
//  be supplied:
//
//     1.
//     void StreamWrite ( RCFile f, RPCClass Object )  {
//       StreamWrite ( f,(PCStream)PCClass );
//     }
//
//     2.
//     PCStream CClassInit ( RPCStream Object )  {
//       return (PCStream)(new CClass(Object));
//     }
//
//     3.
//     void StreamRead ( RCFile f, RPCClass Object )  {
//       StreamRead_ ( f,(PCStream)Object,CClassInit );
//     }
//
//    All these functions are automatically generated by macros
//  DefineStreamFunctions(CClass) -- in the header -- and
//  MakeStreamFunctions(CClass) -- in the implementation body.
//  Then CClass may be streamed in/out using functions #1 and #3.
//    StreamRead will return NULL for Object if it was not
//  in the stream. If Object existed before StreamRead(..) but
//  was not found in the stream, it will be disposed.

void StreamRead_ ( RCFile f, RPCStream Object,
                   InitStreamObject Init )  {
int i;
  f.ReadInt ( &i );
  if (i)  {
    if (!Object)
      Object = Init(Object); //Object = new CStream ( Object );
    Object->read ( f );
  } else  {
    if (Object)  delete Object;
    Object = NULL;
  }
}

void StreamWrite_ ( RCFile f, RPCStream Object )  {
int i;
  if (Object)  {
    i = 1;
    f.WriteInt ( &i );
    Object->write ( f );
  } else  {
    i = 0;
    f.WriteInt ( &i );
  }
}

MakeStreamFunctions(CStream)
