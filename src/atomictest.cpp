#include <atomic>
 #include <new>
 int main(int, char**){
 double a = 1.;
 std::atomic<double> *atom = new(&a) std::atomic<double>;
 double old = *atom;
 while (!atom->compare_exchange_weak(old, old + 1.)) {}
 return 0;}
