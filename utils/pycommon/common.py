#copied from cheng's MPAS-LETKF

import logging
import datetime as dt

def split_timedelta(delt=dt.timedelta(days=0,hours=4,minutes=40,seconds=20)):
    total_seconds = delt.total_seconds()
    days    = int( total_seconds/3600/24 )
    hours   = int( (total_seconds - days*3600*24)/3600 )
    minutes = int( (total_seconds - days*24*3600 - hours*3600)/60 )
    seconds = int( total_seconds - ((days*24 + hours)*60 + minutes)*60 )
    return days, hours, minutes, seconds

def setupEasyLogging(logFilePrefix=None, screenLog=True, fileLog=True):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    date_formatter = '%Y/%m/%d %H:%M:%S'
    #message_formatter  = '[%(asctime)s] %(levelname)s (%(filename)s): %(message)s'
    message_formatter  = '[%(asctime)s] %(levelname)s : %(message)s'
    formatter = logging.Formatter(message_formatter, date_formatter)

    if screenLog:
        hdlr_screen_log = logging.StreamHandler()
        hdlr_screen_log.setFormatter(formatter)
        logger.addHandler(hdlr_screen_log)

    if fileLog:
        if logFilePrefix is None:
            logFilePrefix = "unknown"
        log_file_name = "LOG.{}_{}".format(logFilePrefix, dt.datetime.now().strftime("%Y-%m-%d_%H:%M:%S"))
        hdlr_file_log = logging.FileHandler(log_file_name)
        hdlr_file_log.setFormatter(formatter)
        logger.addHandler(hdlr_file_log)
        #logger.info("log info written into file: {}".format(log_file_name))

    return logger


if __name__ == '__main__':
    #print("{0}_{:2d}:{:2d}:{:2d}".format(split_timedelta()))
    print("{0}_{1:02d}:{2:02d}:{3:02d}".format(*split_timedelta()))
