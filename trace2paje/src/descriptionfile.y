%{
  #include <stdio.h>
  #include <string.h>
  #include "tree.hh"
  #include "container.hh"
  #include "descriptionfileparser.hpp"

  #include <iostream>
  using namespace std;
  using namespace Paje;


  extern "C"
  {
          int yylex(void);
          int yyparse(void);

  }
  int yywrap()
  {
          return 1;
  }
  void yyerror(const char *str)
  {
    fprintf(stderr,"error: %s\n",str);
  }






%}

%union {
  char* string;
  long int_lit;
  void * ptr;
  Paje::Container   * container;
  tree<struct semantic_attribute*>::iterator *attrib_iterator;

}



%token ENDL;
%token<string> IDENTIFIER
%token<string> STRING_LIT
%token<int_lit> NUMBER
%token<ptr> TOK_CONTAINER;
%token<ptr> TOK_RASTRO_HEADER
%token<ptr> TOK_EVENT_TYPE;
%token<ptr> TOK_LINKTYPE ;
%token<ptr> TOK_NAME     ;
%token<ptr> TOK_CREATE   ;
%token<ptr> TOK_PARENT   ;
%token<ptr> TOK_DESTROY  ;
%token<ptr> TOK_CHILDREN ;
%token<ptr> TOK_STATE    ;
%token<ptr> TOK_EVENT    ;
%token<ptr> TOK_ACCEPT   ;
%token<ptr> TOK_IGNORE   ;
%token<ptr> TOK_ALL      ;
%token<ptr> TOK_INCLUDE  ;
%token<ptr> TOK_VALUE    ;
%token<ptr> TOK_TYPE     ;

%type<ptr> code
%type<ptr> header
%type<ptr> definition
%type<ptr> definitions
%type<ptr> container_definition
%type<ptr> toplevel_container_definition
%type<attrib_iterator> container_param
%type<attrib_iterator> container_params
%type<attrib_iterator> name_param
%type<attrib_iterator> create_param
%type<attrib_iterator> destroy_param
%type<string> idf



%%

code:
          header
          definitions
            {$$ = $1}
        ;

header:
        TOK_RASTRO_HEADER '(' idf ',' idf ')'  { 
                  cout <<   "Rastro !" << $3 <<" " << $5 << endl;
                  $$ = NULL;
                  }
        ;


definitions: definition
           | definitions definition
           ;

definition:
          toplevel_container_definition {} 
          ;

toplevel_container_definition:
                             container_definition
                             ;

container_definition:
          TOK_CONTAINER IDENTIFIER '{' container_params '}' {
              attribs_t::iterator *it = $4;

              Paje::Container* c = new Paje::Container(*it);
              print_tree(attributes,attributes->begin(),
              attributes->end());

              reparent_containers(c, *it);

            }
          ;

container_params:
                  container_params container_param {
                    attribs_t::iterator *it;
                    it = new attribs_t::iterator();

                    if ( (*(*$1))->id == ID_NOP ) {
                      attributes->reparent(*$1,*$2);
                      $$ = $1;
                    } else {
                      struct semantic_attribute * nop = new_semantic_attribute();
                      nop->id = ID_NOP;

                      *it = attributes->insert(attributes->begin(),nop);
                      attributes->reparent(*it,*$1);
                      attributes->reparent(*it,*$2);
                      $$ = it;
                    }
                  }
                | container_param { $$ = $1}
                ;

container_param: name_param  
                |create_param 
                |destroy_param 
                | container_definition
                ;

name_param: TOK_NAME STRING_LIT {
                    struct semantic_attribute * attr;
                    attribs_t::iterator *it;

                    it = new attribs_t::iterator();
                    attr = new_semantic_attribute();

                    attr->id = ID_NAME;
                    attr->vals.name = $2;

                    *it = attributes->insert(attributes->begin(),attr);

                    //print_tree(attributes,(*it),(*it).end());
                    //print_tree(attributes,attributes->begin(),attributes->end());

                    $$ = it;
                    cout << "Name: " << $2 <<endl;
                  } ;

          

create_param: TOK_CREATE IDENTIFIER {
                  struct semantic_attribute * attr;
                  attribs_t::iterator *it;

                  it = new attribs_t::iterator();
                  attr = new_semantic_attribute();

                  attr->id = ID_CREATE_EVENT;
                  attr->vals.create_event = $2;

                  *it = attributes->insert(attributes->begin(),attr);

                  //print_tree(attributes,(*it),(*it).end());
                  //print_tree(attributes,attributes->begin(),attributes->end());

                  $$ = it;
                  cout << "Create on event " << $2 <<endl; 
                }
            | TOK_CREATE TOK_PARENT {
                  struct semantic_attribute * attr;
                  attribs_t::iterator *it;

                  it = new attribs_t::iterator();
                  attr = new_semantic_attribute();

                  attr->id = ID_CREATE_PARENT;
                  attr->vals.create_parent = true;

                  *it = attributes->insert(attributes->begin(),attr);

                  //print_tree(attributes,(*it),(*it).end());
                  //print_tree(attributes,attributes->begin(),attributes->end());

                  $$ = it;
                  cout << "Create on parent" <<endl;
                }
            ;


destroy_param: TOK_DESTROY IDENTIFIER {
                  struct semantic_attribute * attr;
                  attribs_t::iterator *it;

                  it = new attribs_t::iterator();
                  attr = new_semantic_attribute();

                  attr->id = ID_DESTROY_EVENT;
                  attr->vals.destroy_event = $2;

                  *it = attributes->insert(attributes->begin(),attr);

                  //print_tree(attributes,(*it),(*it).end());
                  //print_tree(attributes,attributes->begin(),attributes->end());

                  $$ = it;
                  cout << "Destroy on event " << $2 <<endl; 
                }
            | TOK_DESTROY TOK_CHILDREN {
                  struct semantic_attribute * attr;
                  attribs_t::iterator *it;

                  it = new attribs_t::iterator();
                  attr = new_semantic_attribute();

                  attr->id = ID_DESTROY_EVENT;
                  attr->vals.destroy_children = true;

                  *it = attributes->insert(attributes->begin(),attr);

                  //print_tree(attributes,(*it),(*it).end());
                  //print_tree(attributes,attributes->begin(),attributes->end());

                  $$ = it;
                  cout << "Destroy on children" <<endl;
                }
            ;


idf: IDENTIFIER { $$ = $1} ;

