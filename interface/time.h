#ifndef __CUSTOMTIMER__
#define __CUSTOMTIMER__

#include <ctime>

/*! a quick time class for input/output timing info
 */

class Time{
public:
  int hr, min, sec;
  Time(clock_t a=0)
  {
    a/=CLOCKS_PER_SEC;
    hr=a/3600;
    a= a%3600;
    min=a/60;
    sec=a%60;
  }
  //!< default constructor

  friend ostream& operator<<(ostream& os, const Time& obj)
  { 
    if(obj.hr > 0)
      os<<obj.hr<<"hr "<<obj.min<<"m "<<obj.sec<<"s";
    else if(obj.min > 0)
      os<<obj.min<<"m "<<obj.sec<<"s";
    else
      os<<obj.sec<<"s";

    return os;
  }    
  //!< for formatting time

};

#endif
