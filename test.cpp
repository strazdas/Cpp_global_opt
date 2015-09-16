// set_intersection example
#include <iostream>     // std::cout
#include <algorithm>    // std::set_intersection, std::sort
#include <vector>       // std::vector
#include <list>         // std::list

using namespace std;

int main () {
  // int first[] = {5,10,15,20,25};
  // int second[] = {50,40,30,20,10};
  list<int> first_set;
  list<int> second_set;
  first_set.push_back(5);
  first_set.push_back(10);
  first_set.push_back(15);
  first_set.push_back(20);
  first_set.push_back(25);
  second_set.push_back(50);
  second_set.push_back(40);
  second_set.push_back(30);
  second_set.push_back(20);
  second_set.push_back(10);

  vector<int> v(10);                            // 0  0  0  0  0  0  0  0  0  0
  vector<int>::iterator it;

  first_set.sort();     //  5 10 15 20 25
  second_set.sort();    // 10 20 30 40 50
  //// Vector sorting:
  // sort(first_set.begin(), first_set.end());     //  5 10 15 20 25
  // sort(second_set.begin(), second_set.end());   // 10 20 30 40 50

  it=set_intersection (first_set.begin(), first_set.end(), second_set.begin(), second_set.end(), v.begin());
                                                // 10 20 0  0  0  0  0  0  0  0
  v.resize(it-v.begin());                       // 10 20

  cout << "The intersection has " << (v.size()) << " elements:\n";
  for (it=v.begin(); it!=v.end(); ++it)
    cout << ' ' << *it;
  cout << '\n';

  return 0;
}
