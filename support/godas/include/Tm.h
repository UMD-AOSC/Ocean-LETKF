#ifndef _TM_H_
#define _TM_H_

#define RYEAR 1900
#define RDOW 2

typedef struct {
               unsigned short year;
               unsigned short month;
               unsigned short day;
               unsigned short hour;
               unsigned short year0;
               unsigned int yday;
               } YMD;

unsigned int YearDay( YMD );

void CalendarDay( YMD* );

unsigned short DayOfWeek( YMD );

#endif
