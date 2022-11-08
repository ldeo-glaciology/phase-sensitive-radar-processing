# setup_logger.py


def _start_logging_in_file_and_console():

    logger = logging.getLogger("default_logger")
    logger.setLevel(level=logging.WARNING)

    logStreamFormatter = logging.Formatter(
        fmt=f"%(levelname)-8s %(asctime)s \t %(filename)s @function %(funcName)s line %(lineno)s - %(message)s", 
        datefmt="%H:%M:%S"
    )
    consoleHandler = logging.StreamHandler(stream=sys.stdout)
    consoleHandler.setFormatter(logStreamFormatter)
    consoleHandler.setLevel(level=logging.DEBUG)

    logger.addHandler(consoleHandler)


    logFileFormatter = logging.Formatter(
        fmt=f"%(levelname)s %(asctime)s (%(relativeCreated)d) \t %(pathname)s F%(funcName)s L%(lineno)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fileHandler = logging.FileHandler(filename='test4.log')
    fileHandler.setFormatter(logFileFormatter)
    fileHandler.setLevel(level=logging.DEBUG)

    logger.addHandler(fileHandler)

    return logger