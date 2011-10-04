// C++ source code
// File: "/home/kassick/Work/olam/trace2paje/src/paje_functions.hh"
// Created: "Ter, 04 Out 2011 14:17:16 -0300 (kassick)"
// Updated: "Ter, 04 Out 2011 14:17:35 -0300 (kassick)"
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

void init_paje_events();
void paje_header(ostream &out);
void pajeDefineContainerType(string &alias,
                             string &containerType, 
                             string &name,
                             ostream &out);

void pajeDefineStateType(string &alias,
                         string &containerType, 
                         string &name,
                         ostream &out);

void pajeDefineLinkType(string &alias,
                        string &containerType,
                        string &sourceContainerType,
                        string &destContainerType, string &name,
                        ostream &out);

void pajeCreateContainer(double timestamp,
                         string &alias,
                         string &type,
                         string &container, string &name,
                         ostream &out);

void pajeDestroyContainer(double timestamp,
                          string &type, string &container,
                          ostream &out);

void pajeSetState(double timestamp,
                  string &container,
                  string &type, string &value,
                  ostream &out);

void pajePushState(double timestamp,
                   string &container,
                   string &type, string &value,
                   ostream &out);

void pajePopState(double timestamp,
                  string &container, string &type,
                  ostream &out);

void pajeStartLink(double timestamp,
                   string &container,
                   string &type,
                   string &sourceContainer,
                   string &value, string &key,
                   ostream &out);

void pajeEndLink(double timestamp,
                 string &container,
                 string &type,
                 string &endContainer,
                 string &value, string &key,
                 ostream &out);



#endif
