import matplotlib as plt

from configobj import ConfigObj, SimpleVal
from validate import Validator


def check_file_exists(file_name) :
    """Check if the file actually exists or throw exception """
    try:
        with open(file_name) as file:
            pass
    except IOError as e:
        print("Unable to open file: ", file_name) #Does not exist OR no read permissions
        exit(1)


def check_config_validity(configFile, config, configSpecFile, configSpec):
    """Check the validity of the configuration file with respect to its specifications"""
    val = Validator()
    test = config.validate(val, preserve_errors=True)
    if test == True:
        print('Extensive Test of Configuration File ', configFile, ' according to specifications ' , configSpecFile, ' Succeeded.')
    else:
        print('Config file failed specifications test: \n  ', config.validate(val, preserve_errors=True ))
        exit(1)

    val = SimpleVal()
    test = config.validate(val)
    if test == True:
        print('All values present.')
    elif test == False:
        print('No values present!')
    else:
        for entry in test:
            if test[entry] == False:
                print('"{}" missing.'.format(entry))


def check_config_boolean(config):
    if config['ConfigCheck']['checkWorks']:
        print('VALIDATION CHECK SUCCEDED')
    else:
        print('VALIDATION CHECK FAILED')
        exit(1)


def load_config(conf_file, conf_spec_file):

    check_file_exists(conf_file)
    check_file_exists(conf_spec_file)

    ##############################################################################
    # CHECK CONFIG FILE AND LOAD PARAMETERS INTO CONFIG OBJECT
    ##############################################################################
    config_spec = ConfigObj(conf_spec_file, interpolation=False, list_values=False, _inspec=True)
    config = ConfigObj(conf_file, configspec=config_spec)

    # Check validity
    check_config_validity(conf_file, config, conf_spec_file, config_spec)

    check_config_boolean(config)

    return config