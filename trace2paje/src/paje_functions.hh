// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/paje_functions.hh"
// Created: "Ter, 04 Out 2011 14:17:16 -0300 (kassick)"
// Updated: "Sex, 07 Out 2011 15:14:54 -0300 (kassick)"
// $Id$
// Copyright (C) 2011, Rodrigo Virote Kassick <rvkassick@inf.ufrgs.br> 
/*
 * ===========================================================================
 * 
 *       Filename:  paje_functions.hh
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  04-10-2011 14:17:16 BRT
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   (), 
 *        Company:  
 * 
 * ===========================================================================
 */


#ifndef __PAJEFUNCTIONS_H__
#define __PAJEFUNCTIONS_H__


#include <string>
#include <iostream>
#include <sstream>

using namespace std;


#define PAJE_ROOT_CONTAINER  "0"
#define PAJE_TS_PRECISION 6

void init_paje_events();
void paje_header(ostream &out);
void pajeDefineContainerType(const string &alias,
                             const string &containerType, 
                             const string &name,
                             ostream &out);

void pajeDefineStateType(const string &alias,
                         const string &containerType, 
                         const string &name,
                         ostream &out);

void pajeDefineLinkType(const string &alias,
                        const string &containerType,
                        const string &sourceContainerType,
                        const string &destContainerType, const string &name,
                        ostream &out);

void pajeCreateContainer(double timestamp,
                         const string &alias,
                         const string &type,
                         const string &container, const string &name,
                         ostream &out);

void pajeDestroyContainer(double timestamp,
                          const string &type, const string &container,
                          ostream &out);

void pajeSetState(double timestamp,
                  const string &container,
                  const string &type, const string &value,
                  ostream &out);

void pajePushState(double timestamp,
                   const string &container,
                   const string &type, const string &value,
                   ostream &out);

void pajePopState(double timestamp,
                  const string &container, const string &type,
                  ostream &out);

void pajeStartLink(double timestamp,
                   const string &container,
                   const string &type,
                   const string &sourceContainer,
                   const string &value, const string &key,
                   ostream &out);

void pajeEndLink(double timestamp,
                 const string &container,
                 const string &type,
                 const string &endContainer,
                 const string &value, const string &key,
                 ostream &out);



#endif
