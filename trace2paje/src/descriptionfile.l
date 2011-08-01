%{
#include <stdlib.h>
#include <stdio.h>
#include "container.hh"
#include "tree.hpp"

#include "descriptionfileparser.hpp"

#define BUFFER_SIZE 512

extern "C"
{
        int yyline;
        int yylex(void);

}

#define YYDECL int yylex(void)

%}
%%

RASTRO                  {
                          return TOK_RASTRO_HEADER;
                        }
CONTAINER               return TOK_CONTAINER;
EVENTTYPE               return TOK_EVENT_TYPE;
LINKTYPE                return TOK_LINKTYPE ;
NAME                    return TOK_NAME     ;
CREATE                  return TOK_CREATE   ;
PARENT                  return TOK_PARENT   ;
DESTROY                 return TOK_DESTROY  ;
CHILDREN                return TOK_CHILDREN ;
STATE                   return TOK_STATE    ;
EVENT                   return TOK_EVENT    ;
ACCEPT                  return TOK_ACCEPT   ;
IGNORE                  return TOK_IGNORE   ;
INCLUDE                 return TOK_INCLUDE  ;
VALUE                   return TOK_VALUE    ;
TYPE                    return TOK_TYPE     ;

[0-9]+                  {
                          yylval.int_lit = atol(yytext);
                          return NUMBER;
                        }
[\@a-zA-Z][a-zA-Z_0-9]*   {
                          yylval.string = strdup(yytext);
                          return IDENTIFIER;
                        }

"\""([^\n\"\\]*(\\[.\n])*)*"\""  {
                          yylval.string=strdup(yytext);
                          return ( STRING_LIT );
                        }
[(),{}]                 {
                          return yytext[0];
                        }

\n                      { yyline++ ;}
[ \t]+                  /* ignore whitespace */;
.                       {
                          printf("Unexpected token `` %c ''\n",yytext[0]) ;
                        }
%%

#ifdef UNIT_TEST

int main(int argc, char** argv)
{

  while (1) 
  {
    int token = yylex();
    if (token != 0)
      printf("got token %d\n",token);
    else
      break;
  }
}

#endif
