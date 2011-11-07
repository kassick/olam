// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/rastro_generate.h"
// Created: "Seg, 07 Nov 2011 16:30:18 -0200 (kassick)"
// Updated: "Seg, 07 Nov 2011 16:58:04 -0200 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  rastro_generate.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  07-11-2011 16:30:18 BRST
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */


#ifndef  __RASTRO_GENERATE_H__
#define  __RASTRO_GENERATE_H__


#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void generate_rst_functions(char has_csource, const char* csource,
                            char has_cheader, const char* cheader,
                            char has_fmodule, const char* fmodule,
                            FILE * signatures_file);

#ifdef __cplusplus
}
#endif

#endif     /* -----  __RASTRO_GENERATE_H__  ----- */


