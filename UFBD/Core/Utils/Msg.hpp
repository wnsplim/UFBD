#ifndef Msg_H
#define Msg_H

#include <string>



namespace Msg {

   void   error(std::string s);
   void   warning(std::string s);
   void   setDeferWarnings(bool b);
   bool   hasDeferredWarnings(void);
   void   flushWarnings(void);
}

#endif
