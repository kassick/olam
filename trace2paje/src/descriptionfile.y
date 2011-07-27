%{
  #include <stdio.h>
  #include <string.h>
  #include "tree.hh"
  #include "container.hh"
  #include "descriptionfileparser.hpp"

  #include <iostream>
  using namespace std;

  typedef tree<Paje::Container*> hierarchy_t;

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


  hierarchy_t * toplevel_hierarchy;




  void init_desc_parser() {
    toplevel_hierarchy = new tree<Paje::Container*>();
  }

%}

%union {
  char* string;
  long int_lit;
  void * ptr;
}



%token<string> IDENTIFIER
%token<string> STRING_LIT
%token<int_lit> NUMBER
%token<ptr> TOK_RASTRO_HEADER
%token<ptr> TOK_CONTAINER;
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
%type<ptr> container_param
%type<ptr> container_params
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
          container_definition { 
              toplevel_hierarchy->insert_subtree(toplevel_hierarchy->begin(), $1)
            }
          ;

container_definition:
          TOK_CONTAINER IDENTIFIER '{' container_params '}'  { 
              cout << "Creating container " << $2 <<endl;
              Paje::Container * c = new Paje::Container();
              c->typeName = $2;

              $$ = new tree<Paje::Container*>();
              $$->set_head(c);

            }
          ;



container_params: container_param
                | container_params container_param
                ;

container_param: name_param
                |create_param
                |destroy_param
                |container_definition
                ;

name_param: TOK_NAME STRING_LIT { cout << "Name: " << $2 <<endl; } ;

create_param: TOK_CREATE IDENTIFIER { 
                  cout << "Create on event " << $2 <<endl; 
                }
            | TOK_CREATE TOK_PARENT {
                  cout << "Create on parent" <<endl;
                }
            ;


destroy_param: TOK_DESTROY IDENTIFIER { 
                  cout << "Destroy on event " << $2 <<endl; 
                }
            | TOK_DESTROY TOK_CHILDREN {
                  cout << "Destroy on children" <<endl;
                }
            ;


idf: IDENTIFIER { $$ = $1} ;

