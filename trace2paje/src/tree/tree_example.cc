/* 

	Cadabra: a field-theory motivated computer algebra system.
	Copyright (C) 2001-2009  Kasper Peeters <kasper.peeters@aei.mpg.de>

   This program is free software: you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation, either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
*/

#include <algorithm>
#include <string>
#include <iostream>
#include "tree.hh"

using namespace std;

void print_tree(const tree < std::string > &tr,
                tree < std::string >::iterator it,
                tree < std::string >::iterator end)
{
  if (!tr.is_valid(it))
  {
    cout << "invalid\n" ;
    return;
  }
  int rootdepth = tr.depth(it);
  std::cout << "-----" << std::endl;
  while (it != end) {
    for (int i = 0; i < tr.depth(it) - rootdepth; ++i)
      std::cout << "  ";
    std::cout << (*it) << std::endl << std::flush;
    ++it;
  }
  std::cout << "-----" << std::endl;
}



int main(int, char **)
{
  tree < string > tr;
  tree < string >::iterator top, one, two, loc, banana, first,last;

  top = tr.begin();
  one = tr.insert(top, "mytop");
  first = tr.insert(top, "mytop2");
  first = tr.insert(top, "mytop3");
  banana = tr.insert(one, "111");
  tr.insert(one, "222");
  tr.insert(one, "last?");
  two = tr.append_child(one,"child1?");
  two = tr.append_child(banana,"child2?");
  print_tree(tr, tr.begin(), tr.end());

  cout <<"here" <<endl <<flush;
  first = tr.begin_fixed(top.begin(),1);
  cout <<"here1" <<endl <<flush;

  while (tr.is_valid(first))
  {
    cout << *first <<endl;
    ++first;
  }
  cout <<"here2" <<endl <<flush;
  
    
  //tree < string >::fixed_depth_iterator f = tr.begin(top);
  //print_tree(tr,sib,tr.end(top));
  /*
     one = tr.insert(top, "one");
     two = tr.append_child(one, "two");
     tr.append_child(two, "apple");
     banana = tr.append_child(two, "banana");
     tr.append_child(banana, "cherry");
     tr.append_child(two, "peach");
     tr.append_child(one, "three"); 

  loc = find(tr.begin(), tr.end(), "two");
  loc = tr.begin();
  if (loc != tr.end()) {
    tree < string >::sibling_iterator sib = tr.begin(loc);
    while (sib != tr.end(loc)) {
      cout << (*sib) << endl;
      ++sib;
    }
    cout << endl;
    tree < string >::iterator sib2 = tr.begin(loc);
    tree < string >::iterator end2 = tr.end(loc);
    while (sib2 != end2) {
      for (int i = 0; i < tr.depth(sib2) - 2; ++i)
        cout << " ";
      cout << (*sib2) << endl;
      ++sib2;
    }
  }
  */
}
