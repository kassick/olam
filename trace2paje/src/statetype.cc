// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/statetype.cc"
// Created: "Qui, 06 Out 2011 15:13:16 -0300 (kassick)"
// Updated: "Qui, 06 Out 2011 16:14:50 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 *
 *       Filename:  statetype.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06-10-2011 15:13:16 BRT
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:   (), 
 *        Company:  
 *
 * ==========================================================================
 */


#include "statetype.hh"
#include "container.hh"
#include "paje_functions.hh"

#include <string>
#include <iostream>


void Paje::StateType::do_header(ostream &out)
{
  if (!container)
    cerr << "no container here!!!" << endl;
  pajeDefineStateType(typeName, container->typeName, typeName, out);
}


Paje::StateType::StateType(const string & typeName): Paje::BaseEventType(typeName) {}

Paje::StateType::StateType(const string & typeName,Paje::Container * c):  Paje::BaseEventType(typeName,c) {}


