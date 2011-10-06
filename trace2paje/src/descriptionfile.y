%{
  #include <stdio.h>
  #include <string.h>
  #include "paje.hh"
  #include "tree.hpp"
  #include "container.hh"
  #include "semantics.hh"
  #include "descriptionfileparser.hpp"
  #include "symbols.hh"

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


attribs_t * create_nop_attr(attribs_t * n)
{
  SemanticAttribute * attr = new_semantic_attribute();
  attr->id = ID_NOP;

  attribs_t * t = new attribs_t(attr);
  if (n)
    t->addChild(n);

  return t;
}



%}

%union {
  char* string;
  long long int_lit;
  void * ptr;
  Paje::Container   * container;
  attribs_t * attr_node;
  //tree<struct semantic_attribute*>::iterator *attrib_iterator;

}



%token ENDL;
%token<string> IDENTIFIER ;
%token<string> STRING_LIT ;
%token<int_lit> NUMBER   ;
%token<string> RST_TYPE  ;
%token<ptr> TOK_CONTAINER;
%token<ptr> TOK_RASTRO_HEADER ;
%token<ptr> TOK_EVENT_TYPE;
%token<ptr> TOK_STATE_TYPE;
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
%token<ptr> TOK_ID       ;
%token<ptr> TOK_START    ;
%token<ptr> TOK_END      ;
%token<ptr> TOK_SOURCE   ;
%token<ptr> TOK_DEST     ;
%token<ptr> TOK_LINK     ;
%token<ptr> TOK_KEY      ;


%type<ptr> code
%type<ptr> header
%type<ptr> definition
%type<ptr> definitions
%type<attr_node> container_definition
%type<attr_node> container_param
%type<attr_node> container_params
%type<attr_node> name_param
%type<attr_node> create_param
%type<attr_node> destroy_param
%type<string> idf
%type<attr_node> accept_item;
%type<attr_node> accept_list;
%type<attr_node> opt_accept_list;
%type<attr_node> accept_param;
%type<attr_node> eventtype_param;
%type<attr_node> statetype_param;
%type<ptr> include_statement;

%type<attr_node> idf_list;
%type<attr_node> idf_idlist;
%type<attr_node> idf_list_item;

%type<attr_node> event_statement;
%type<attr_node> event_params;
%type<attr_node> event_param;
%type<attr_node> state_param;
%type<attr_node> state_params;
%type<attr_node> state_statement;
%type<attr_node> link_param;
%type<attr_node> link_params;
%type<attr_node> link_statement;
%type<attr_node> event_id_statement;
%type<attr_node> linktype_param;



%%

code:
          header
          definitions
            {$$ = $1}
        ;

header:
        TOK_RASTRO_HEADER '(' idf ',' idf ')'  { 
                  //cout <<   "Rastro !" << $3 <<" " << $5 << endl;
                  Paje::idf1_name = $3;
                  Paje::idf2_name = $5;
                  $$ = NULL;
                  }
        | { // In the case of files included, we do not expect them to have
            // RASTRO definitions
            $$ = NULL;
          }
        ;


definitions: definition
           | definitions definition
           ;

definition:
           container_definition {
              early_parse_tree->addChild($1);
           }
          |include_statement {}
          |event_statement {
            early_parse_tree->addChild($1);
          }
          |state_statement {
            early_parse_tree->addChild($1);
          }
          |link_statement {
            early_parse_tree->addChild($1);
          }
          |event_id_statement { 
              late_parse_tree->addChild ($1);
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
                    /*$$ = create_nop_attr(NULL);
                    $$->addChild($1);
                    $$->addChild($2);*/
                    //$1->addChild(create_nop_attr($2));
                    $1->addChild($2);
                  }
                | container_param { $$ = create_nop_attr($1);}
                ;

container_param: name_param
                | create_param
                | destroy_param
                | accept_param
                | eventtype_param
                | linktype_param
                | statetype_param
                | container_definition {
                  $$ = create_nop_attr($1);
                  }
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
            | TOK_DESTROY TOK_PARENT {
                  attr = new_semantic_attribute();
                  attr->id = ID_DESTROY_WITH_PARENT;
                  attr->vals.destroy_with_parent = true;

                  n = new attribs_t(attr);
                  $$ = n;

                  //cout << "Destroy on children" <<endl;
                }
            ;



statetype_param: TOK_STATE_TYPE IDENTIFIER '{' accept_list '}' {
                  attr = new_semantic_attribute();
                  attr->id = ID_STATE_TYPE_DEF;
                  attr->vals.name = $2;

                  n = new attribs_t(attr);
                  n->addChild($4);
                  $$ = n;
                  
               }
               | TOK_STATE_TYPE IDENTIFIER {
                  attr = new_semantic_attribute();
                  attr->id = ID_STATE_TYPE_DEF;
                  attr->vals.name = $2;

                  n = new attribs_t(attr);
                  $$ = n;
                  
               }


eventtype_param: TOK_EVENT_TYPE IDENTIFIER '{' accept_list '}' {
                  attr = new_semantic_attribute();
                  attr->id = ID_EVENT_TYPE_DEF;
                  attr->vals.name = $2;

                  n = new attribs_t(attr);
                  n->addChild($4);
                  $$ = n;
                  
               }
               | TOK_EVENT_TYPE IDENTIFIER {
                  attr = new_semantic_attribute();
                  attr->id = ID_EVENT_TYPE_DEF;
                  attr->vals.name = $2;

                  n = new attribs_t(attr);
                  $$ = n;
                  
               }

linktype_param: TOK_LINKTYPE IDENTIFIER TOK_SOURCE IDENTIFIER TOK_DEST
                            IDENTIFIER opt_accept_list  {
                  // link type node with type name
                  attr = new_semantic_attribute();
                  attr->id = ID_LINK_TYPE;
                  attr->vals.name = $2;

                  n = new attribs_t(attr);

                  // source node
                  attr = new_semantic_attribute();
                  attr->id = ID_LINK_SOURCE;
                  attr->vals.name = $4;
                  n->addChild(new attribs_t(attr));
                  
                  // dest node
                  attr = new_semantic_attribute();
                  attr->id = ID_LINK_DEST;
                  attr->vals.name = $6;
                  n->addChild(new attribs_t(attr));
                  
                  if ($7) {
                    // add the accept list
                    n->addChild($7);
                  }

                  $$ = n;
                }

opt_accept_list: '{' accept_list '}' {
               $$ = $2;
               }
               | {$$ = NULL};


accept_param: TOK_ACCEPT '{' accept_list '}' {

              //cout << "Accept! " <<endl;
              //print_tree($3);

              $$ = $3;
              }
            | TOK_ACCEPT '{' '}' {

              //cout << "Accept! " <<endl;
              //print_tree($3);

              $$ = NULL;
            }




accept_list: accept_item {$$ = $1;}
           | accept_list accept_item {
                    $1->addChild($2);
                    $$ = $1;
                  }
          ;

accept_item: IDENTIFIER {
                    attr = new_semantic_attribute();
                    attr->id = ID_ACCEPT_LIST;
                    attr->vals.identifier_name = $1;

                    $$ = new attribs_t(attr);
                  }
          ;




idf_list: no_commas idf_idlist no_commas { $$ = $2; }
        ;

idf_idlist: idf_list_item {$$ = $1;}
          | idf_idlist one_comma idf_list_item {
                  $1->addChild($3);
                  $$ = $1;
              }
          ;
one_comma: commas
         |                 {}
         ;
no_commas: commas          {}
         |
         ;
commas   : commas ','      {}
         | ','
         ;

idf_list_item: IDENTIFIER {
                    attr = new_semantic_attribute();
                    attr->id = ID_IDF;
                    attr->vals.identifier_name = $1;

                    $$ = new attribs_t(attr);
                  }
          ;



idf: IDENTIFIER { $$ = $1} ;




include_statement: TOK_INCLUDE STRING_LIT {
                  $$ = NULL;

                  files_to_parse.push($2);
                }
                
event_statement: TOK_EVENT IDENTIFIER '{' event_params '}' {
                  attr = new_semantic_attribute();
                  attr->id = ID_EVENT;
                  attr->vals.name = $2;

                  n = new attribs_t(attr);
                  n->addChild($4);
                  $$ = n;
               }

event_params: event_params event_param {
                    $1->addChild($2);
                    $$ = $1;
              }
            | event_param {
              $$ = create_nop_attr($1);
            }
            ;


event_param: TOK_TYPE IDENTIFIER {
                  attr = new_semantic_attribute();
                  attr->id = ID_EVENT_TYPE;
                  attr->vals.name = $2;

                  n = new attribs_t(attr);
                  $$ = n;

               }
          | RST_TYPE IDENTIFIER {
                  attr = new_semantic_attribute();
                  attr->id = ID_RASTRO_TYPE;
                  attr->vals.name = $1;

                  n = new attribs_t(attr);

                  SemanticAttribute * attr2 = new_semantic_attribute();
                  attr2->id = ID_IDF;
                  attr2->vals.name = $2;

                  attribs_t *n2 = new attribs_t(attr2);

                  n->addChild(n2);

                  $$ = n;
                }
          | RST_TYPE '{' idf_list '}' {
                  attr = new_semantic_attribute();
                  attr->id = ID_RASTRO_TYPE;
                  attr->vals.name = $1;

                  n = new attribs_t(attr);
                  n->addChild($3);

                  $$ = n;
                }
          | TOK_VALUE STRING_LIT {

                  attr = new_semantic_attribute();
                  attr->id = ID_FORMAT_NAME;
                  attr->vals.name = $2;
                  
                  n = new attribs_t(attr);
                  $$ = n;
            }
          ;

state_statement: TOK_STATE IDENTIFIER '{' state_params '}' {
            attr = new_semantic_attribute();
            attr->id = ID_STATE;
            attr->vals.name = $2;

            n = new attribs_t(attr);
            n->addChild($4);
            $$ = n;
            /*
            // create a new state instance
            // call fill_from_attr
            Paje::State *state ;
            if (!(event_names->count($2)))
            {
              string name($2);
              state = new Paje::State(name);
              (*event_names)[$2] = state;
            } else {
              state = (Paje::State * )(*event_names)[$2];
            }
            state->fill_from_attr($4);

            $$ = NULL; */
          }

state_params: state_params state_param    {
                    $1->addChild($2);
                    $$ = $1;
              }
            | state_param {$$ = create_nop_attr($1);}


state_param: TOK_TYPE IDENTIFIER {
                  attr = new_semantic_attribute();
                  attr->id = ID_EVENT_TYPE;
                  attr->vals.name = $2;

                  n = new attribs_t(attr);
                  $$ = n;

               }
          | RST_TYPE IDENTIFIER {
                  attr = new_semantic_attribute();
                  attr->id = ID_RASTRO_TYPE;
                  attr->vals.name = $1;

                  n = new attribs_t(attr);

                  SemanticAttribute * attr2 = new_semantic_attribute();
                  attr2->id = ID_IDF;
                  attr2->vals.name = $2;

                  attribs_t *n2 = new attribs_t(attr2);

                  n->addChild(n2);

                  $$ = n;
                }
          | RST_TYPE '{' idf_list '}' {
                  attr = new_semantic_attribute();
                  attr->id = ID_RASTRO_TYPE;
                  attr->vals.name = $1;

                  n = new attribs_t(attr);
                  n->addChild($3);

                  $$ = n;
                }
          | TOK_VALUE STRING_LIT {
                  // remove ""

                  attr = new_semantic_attribute();
                  attr->id = ID_FORMAT_NAME;
                  attr->vals.name = $2;
                  
                  n = new attribs_t(attr);
                  $$ = n;
            }
          ;



link_statement: TOK_LINK IDENTIFIER '{' link_params '}' {
            // create a new state instance
            // call fill_from_attr

            attr = new_semantic_attribute();
            attr->id = ID_LINK;
            attr->vals.name = $2;

            n = new attribs_t(attr);
            n->addChild($4);
            $$ = n;
          }

link_params: link_params link_param    {
                    $1->addChild($2);
                    $$ = $1;
              }
            | link_param {$$ = create_nop_attr($1);}


link_param: TOK_TYPE IDENTIFIER {
                  attr = new_semantic_attribute();
                  attr->id = ID_EVENT_TYPE;
                  attr->vals.name = $2;

                  n = new attribs_t(attr);
                  $$ = n;

               }
          | RST_TYPE IDENTIFIER {
                  attr = new_semantic_attribute();
                  attr->id = ID_RASTRO_TYPE;
                  attr->vals.name = $1;

                  n = new attribs_t(attr);

                  SemanticAttribute * attr2 = new_semantic_attribute();
                  attr2->id = ID_IDF;
                  attr2->vals.name = $2;

                  attribs_t *n2 = new attribs_t(attr2);

                  n->addChild(n2);

                  $$ = n;
                }
          | RST_TYPE '{' idf_list '}' {
                  attr = new_semantic_attribute();
                  attr->id = ID_RASTRO_TYPE;
                  attr->vals.name = $1;

                  n = new attribs_t(attr);
                  n->addChild($3);

                  $$ = n;
                }
          | TOK_VALUE STRING_LIT {

                  attr = new_semantic_attribute();
                  attr->id = ID_FORMAT_NAME;
                  attr->vals.name = $2;
                  
                  n = new attribs_t(attr);
                  $$ = n;
            }
          | TOK_KEY STRING_LIT {

                  attr = new_semantic_attribute();
                  attr->id = ID_KEY_FORMAT;
                  attr->vals.name = $2;
                  
                  n = new attribs_t(attr);
                  $$ = n;
            }
          ;


event_id_statement:
              TOK_ID IDENTIFIER NUMBER {
                // simple event
                attr = new_semantic_attribute();
                attr->id = ID_EVENT_START;
                attr->vals.name = $2;

                n = new attribs_t(attr);
                
                attr = new_semantic_attribute();
                attr->id = ID_EVENT_ID;
                attr->vals.event_id = $3;

                attribs_t * n1 = new attribs_t(attr);
                n->addChild(n1);
                $$ = n;
              } 
            | TOK_ID IDENTIFIER TOK_START NUMBER {
                // Start of state 
                attr = new_semantic_attribute();
                attr->id = ID_STATE_START;
                attr->vals.name = $2;

                n = new attribs_t(attr);
                
                attr = new_semantic_attribute();
                attr->id = ID_EVENT_ID;
                attr->vals.event_id = $4;

                attribs_t * n1 = new attribs_t(attr);
                n->addChild(n1);
                $$ = n;
              }
            | TOK_ID IDENTIFIER TOK_END NUMBER {
                // Start of state 
                attr = new_semantic_attribute();
                attr->id = ID_STATE_END;
                attr->vals.name = $2;

                n = new attribs_t(attr);

                attr = new_semantic_attribute();
                attr->id = ID_EVENT_ID;
                attr->vals.event_id = $4;

                attribs_t * n1 = new attribs_t(attr);
                n->addChild(n1);
                $$ = n;
              }
            | TOK_ID IDENTIFIER NUMBER one_comma NUMBER {
                // Create nop with two childs:
                //  nop-+ 
                //      |
                //      +----START-+
                //      |          |
                //      |          +-ID
                //      +-END-+
                //            |
                //            +-ID


                //      toplevel: nop
                attr = new_semantic_attribute();
                attr->id = ID_NOP;

                n = new attribs_t(attr);
                
                // create the start event
                attr = new_semantic_attribute();
                attr->id = ID_STATE_START;
                attr->vals.name = $2;

                attribs_t *n1 = new attribs_t(attr);

                attr = new_semantic_attribute();
                attr->id = ID_EVENT_ID;
                attr->vals.event_id = $3;

                attribs_t * n2 = new attribs_t(attr);

                n1->addChild(n2);
                n->addChild(n1);
                
                // create the end event
                attr = new_semantic_attribute();
                attr->id = ID_STATE_END;
                attr->vals.name = strdup($2); // otherwise we'll have double free

                n1 = new attribs_t(attr);

                attr = new_semantic_attribute();
                attr->id = ID_EVENT_ID;
                attr->vals.event_id = $5;

                n2 = new attribs_t(attr);

                n1->addChild(n2);
                n->addChild(n1);

                $$ = n;
              }
            ;

