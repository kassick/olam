// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro_loop.hh"
// Created: "Ter, 27 Set 2011 15:33:14 -0300 (kassick)"
// Updated: "Qua, 05 Out 2011 18:43:27 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  rastro_loop.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  27-09-2011 15:33:14 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */


#ifndef __RASTRO_LOOP_HH__
#define __RASTRO_LOOP_HH__

#include "semantics.hh"
#include "paje.hh"
#include "rastro_helper.hh"
#include <ostream>
#include <list>

extern "C" {
#include <rastro.h>
}

#define _RST_BUF_SIZE 100000

double  rastro_loop_events(list<string> &files_to_open, ostream &out) ;


#endif
