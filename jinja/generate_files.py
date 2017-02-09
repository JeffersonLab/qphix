#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

import argparse
import configparser

import jinja2


def main():
    isas = configparser.ConfigParser()
    isas.read('isa.ini')

    kernels = [
        #'clov_dslash',
        #'dslash',
        ('tmf_dslash', 'tmf_dslash'),
        ('tmf_clov_dslash', 'tmf_clov_%(fptype)s_dslash'),
    ]

    # Setting up Jinja
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader('..')
    )
    complete_specialization = env.get_template('jinja/complete_specialization.j2.h')

    for kernel, kernel_pattern in kernels:
        for isa in isas.sections():

            defines = []
            for fptype in isas[isa]:
                val = isas[isa][fptype]

                veclen_str, soalen_str = [x.strip() for x in val.split(';')]
                veclen = int(veclen_str)
                soalens = [int(x.strip()) for x in soalen_str.split()]

                defines.append((fptype, veclen, soalens))

            print(defines)

            rendered = complete_specialization.render(
                isa=isa,
                kernel_base=kernel,
                kernel_pattern=kernel_pattern,
                defines=defines,
            )

            filename = '{}_{}_complete_specialization.h'.format(kernel, isa)
            with open(filename, 'w') as f:
                f.write(rendered)


def _parse_args():
    '''
    Parses the command line arguments.

    :return: Namespace with arguments.
    :rtype: Namespace
    '''
    parser = argparse.ArgumentParser(description='')
    options = parser.parse_args()

    return options


if __name__ == "__main__":
    main()
