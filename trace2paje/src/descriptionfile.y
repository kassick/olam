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
  tree<Paje::Container*>::iterator hierarchy_iterator;
  Paje::Container* container;

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
%type<hierarchy_iterator> container_definition
%type<container> container_param
%type<container> container_params
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
          toplevel_container_definition 
          ;

toplevel_container_definition:
                             {$1 = tr->begin(); } container_definition 
                             ;

container_definition:
          TOK_CONTAINER IDENTIFIER {
                //bottom up approach here; send the iterators
                $3 = new Paje::Container();
                $3->typeName = $2;
                $3->where = toplevel_hierarchy->append_child($$,c);
            } '{' container_params '}'
          ;

container_params: {$1 = $$} container_param
                | {$1 = $2 = $$;} container_params container_param
                ;

container_param: name_param
                |create_param
                |destroy_param
                |{$1 = $$->where } container_definition
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

