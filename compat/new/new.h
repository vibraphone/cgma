#include <new>

// According to section 18.4 of the 2005-10-19 working draft
// of the C++ Standard (Doc no: N1905=05-0165), these are the 
// symbols defined in the <new> header. 

using std::bad_alloc;
using std::nothrow_t;
using std::nothrow;
using std::new_handler;
using std::set_new_handler;
