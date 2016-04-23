#ifndef _TMX_H_
#define _TMX_H_

#define RYEAR 1900
#define RDOW 2

typedef struct {
               int year;
               int month;
               int day;
               int hour;
               int year0;
               unsigned int yday;
               } YMD;

unsigned int YearDay( YMD );

void CalendarDay( YMD* );

int DayOfWeek( YMD );

#endif
