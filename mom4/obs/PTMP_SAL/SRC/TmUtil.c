/*   Time utilities                                 */

#include "Tm.h"

static unsigned short dpm[12] = {
                   31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
unsigned short i, yr, mn;
unsigned int dy, yd;

unsigned int YearDay(YMD date)
{
  yd = 0;
  for (i = 0; i < date.month-1; i++)
    yd += dpm[i];
  yd += date.day;
  if (date.month > 2 && !(date.year%4) && ((date.year%100) || !(date.year%400)))
    yd++;

  for (i = date.year0; i < date.year; i++)
    if (i%4)
      yd += 365;
    else
      if ((i%100) || !(i%400))
        yd += 366;
      else
        yd += 365;

  return(yd);
}

/* ======================================================== */

void CalendarDay(YMD *date)
{
  yr = date->year0;
  dy = date->yday;

  if (yr%4)
    yd = 365;
  else
    if ((yr%100) || !(yr%400))
      yd = 366;
    else
      yd = 365;

  while (yd < dy) {
    dy -= yd;
    yr++;
    if (yr%4)
      yd = 365;
    else
      if ((yr%100) || !(yr%400))
        yd = 366;
      else
        yd = 365;
  }

  if (yr%4)
    dpm[1] = 28;
  else
    if ((yr%100) || !(yr%400))
      dpm[1] = 29;
    else
      dpm[1] = 28;
  mn = 0;
  while (dpm[mn] < dy) {
    dy -= dpm[mn];
    mn++;
  }
  dpm[1] = 28;
  date->year = yr;
  date->month = mn + 1;
  date->day = dy;
}

/* ======================================================== */

unsigned short DayOfWeek(YMD date)
{
  date.year0 = RYEAR;
  return((unsigned short)((YearDay(date)+RDOW-2)%7 + 1));
}
