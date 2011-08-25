%{
  #include <stdio.h>
  #include <string.h>
  #include "paje.hh"
  #include "tree.hpp"
  #include "container.hh"
  #include "semantics.hh"
  #include "descriptionfileparser.hpp"

  #include <iostream>
  using namespace std;
  using namespace Paje;


  extern "C"
  {
          extern int yyline;
          int yylex(void);
          int yyparse(void);

  }
  int yywrap()
  {
          return 1;
  }
  void yyerror(const char *str)
  {
    fprintf(stderr,"error in line %d: %s\n",yyline,str);
  }


SemanticAttribute * attr;
attribs_t *n;



%}

%union {
  char* string;
  long int_lit;
  void * ptr;
  Paje::Container   * container;
  attribs_t * attr_node;
  //tree<struct semantic_attribute*>::iterator *attrib_iterator;

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
%type<attr_node> container_definition
%type<container> toplevel_container_definition
%type<attr_node> container_param
%type<attr_node> container_params
%type<attr_node> name_param
%type<attr_node> create_param
%type<attr_node> destroy_param
%type<string> idf
%type<attr_node> accept_item;
%type<attr_node> accept_list;
%type<attr_node> accept_param;



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
                             container_definition {
                                hierarchy_t * top = attr_to_container_hierarchy($1,toplevel_hierarchy);
                                
                                //cout << "toplevel hierarchy:" <<endl;

                                //print_tree(top);

                                //cout <<"----?" <<endl;

                                //toplevel_hierarchy->addChild(top);
                              
                                //cout << "Defined a toplevel container:" << $1->getVal()->vals.container_type <<endl;


                                //free_subtree($1,true);
                                $$ = NULL;
                              }
                             ;

container_definition:
          TOK_CONTAINER IDENTIFIER '{' container_params '}' {
            //cout << "Container " << $2 << " has attributes : " << endl;
            //print_tree<SemanticAttribute *>($4);
            //cout << "----" <<endl;



            //Container * c = new Container($2,$4);
            attr = new SemanticAttribute();
            attr->id = ID_CONTAINER;
            //attr->vals.container = c;
            attr->vals.container_type = $2;
            $$ = new attribs_t(attr);
            $$->addChild($4);
            
            //cout << "Finished container_definition" <<endl;
          }
          ;

container_params:
                  container_params container_param {
                    $1->addChild($2);
                  }
                | container_param { $$ = $1}
                ;

container_param: name_param
                | create_param
                | destroy_param
                | accept_param
                | container_definition
                ;

name_param: TOK_NAME STRING_LIT {
                    attr = new_semantic_attribute();
                    attr->id = ID_NAME;
                    attr->vals.name = $2;

                    n = new attribs_t(attr);

                    $$ = n;
                    //cout << "Name: " << $2 <<endl;
                  } ;

          

create_param: TOK_CREATE IDENTIFIER {
                  attr = new_semantic_attribute();
                  attr->id = ID_CREATE_EVENT;
                  attr->vals.create_event = $2;

                  n = new attribs_t(attr);
                  $$ = n;

                  //cout << "Create on event " << $2 <<endl; 
                }
            | TOK_CREATE TOK_PARENT {
                  attr = new_semantic_attribute();
                  attr->id = ID_CREATE_PARENT;
                  attr->vals.create_parent = true;

                  n = new attribs_t(attr);
                  $$ = n;
                  //cout << "Create on parent" <<endl;
                }
            ;


destroy_param: TOK_DESTROY IDENTIFIER {
                  attr = new_semantic_attribute();
                  attr->id = ID_DESTROY_EVENT;
                  attr->vals.destroy_event = $2;

                  n = new attribs_t(attr);
                  $$ = n;

                  //cout << "Destroy on event " << $2 <<endl; 
                }
            | TOK_DESTROY TOK_CHILDREN {
                  attr = new_semantic_attribute();
                  attr->id = ID_DESTROY_CHILDREN;
                  attr->vals.destroy_children = true;

                  n = new attribs_t(attr);
                  $$ = n;

                  //cout << "Destroy on children" <<endl;
                }
            ;


accept_param: TOK_ACCEPT '{' accept_list '}' {

              //cout << "Accept! " <<endl;
              //print_tree($3);

              $$ = $3;
            }
            ;

accept_list: accept_list accept_item {
                    $1->addChild($2);
                    $$ = $1;
                  }
               | accept_item {$$ = $1;}
               ;

accept_item: IDENTIFIER {
                    attr = new_semantic_attribute();
                    attr->id = ID_ACCEPT_LIST;
                    attr->vals.identifier_name = $1;

                    $$ = new attribs_t(attr);
                  }
          ;


idf: IDENTIFIER { $$ = $1} ;

