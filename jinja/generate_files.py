#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

import argparse
import getpass
import json
import os
import itertools
import socket

import jinja2


def get_kernel_files_for_isa(kernel_pattern, isa, fptypes):
    generated_files = os.listdir(os.path.join('..', 'generated', isa,
                                              'generated'))
    files = []
    for fptype in fptypes:
        prefix = kernel_pattern % {'fptype_underscore': fptype + '_'}
        added = list(filter(lambda x: x.startswith(prefix), generated_files))
        files += added
    return sorted(set(files))


def main():
    options = _parse_args()

    generated_warning = 'This file has been automatically generated. Do not change it manually, rather look for the template in qphix-codegen.'

    with open('isa.js') as f:
        isas = json.load(f)

    with open('kernels.js') as f:
        kernel_patterns = json.load(f)


    # Setting up Jinja
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader('..')
    )
    complete_specialization = env.get_template('jinja/complete_specialization.h.j2')
    makefile_am = env.get_template('jinja/Makefile.am.j2')
    kernel_generated_h = env.get_template('jinja/kernel_generated.h.j2')

    for kernel_pattern in kernel_patterns:
        kernel = kernel_pattern % {'fptype_underscore': ''}

        rendered = kernel_generated_h.render(
            generated_warning=generated_warning,
            kernel=kernel,
            isas=sorted(isas.keys()),
        )
        filename = '../generated/{}_generated.h'.format(kernel)
        with open(filename, 'w') as f:
            f.write(rendered)

    for isa, isa_data in sorted(isas.items()):
        if len(options.isa) > 0 and not isa in options.isa:
            continue

        if not os.path.isdir(os.path.join('..', 'generated', isa)):
            print('Code for ISA `{}` is not generated. Skipping.'.format(isa))
            continue

        os.makedirs(os.path.join('..', 'generated', isa, 'include'), exist_ok=True)
        os.makedirs(os.path.join('..', 'generated', isa, 'lib'), exist_ok=True)

        # Generate a `Makefile.am`.
        generated_files = os.listdir(os.path.join('..', 'generated', isa, 'generated'))
        generated_for_prefix = {}
        for kernel_pattern in kernel_patterns:
            prefix = kernel_pattern % {'fptype_underscore': ''}
            generated_for_prefix[prefix] = get_kernel_files_for_isa(kernel_pattern, isa, isa_data['fptypes'].keys())

        rendered = makefile_am.render(
            generated_warning=generated_warning,
            isa=isa,
            extra_includes=[
                os.path.relpath(extra_include, os.path.join('qphix', isa))
                for extra_include in isa_data['extra_includes_local']],
            generated_for_prefix=generated_for_prefix,
        )
        filename = os.path.join('../generated', isa, 'Makefile.am')
        with open(filename, 'w') as f:
            f.write(rendered)

        for kernel_pattern in kernel_patterns:
            kernel = kernel_pattern % {'fptype_underscore': ''}

            filename_decl = os.path.join('..', 'generated', isa, 'include', '{}_{}_decl.h'.format(kernel, isa))
            template_decl = env.get_template('jinja/{}_general.h.j2'.format(kernel))

            rendered = template_decl.render()
            with open(filename_decl, 'w') as f:
                f.write(rendered)

            for fptype, fptype_data in sorted(isa_data['fptypes'].items()):
                veclen = fptype_data['veclen']
                for soalen in fptype_data['soalens']:
                    print('Working on kernel `{}` for ISA `{}` …'.format(kernel, isa))
                    for compress12 in ['true', 'false']:
                        for tbc in itertools.product(*([(True, False)] * 4)):
                            tbc_ts = ''.join(['t' if x else 's' for x in tbc])
                            tbc_true_false = ''.join(['true' if x else 'false' for x in tbc])

                            # Generate the complete specialization.
                            rendered = complete_specialization.render(
                                FPTYPE=fptype,
                                VEC=veclen,
                                SOA=soalen,
                                COMPRESS12=compress12,
                                generated_warning=generated_warning,
                                ISA=isa,
                                tbc_ts=tbc_ts,
                                tbc_true_false=tbc_true_false,
                                kernel_base=kernel,
                                kernel_pattern=kernel_pattern,
                                extra_includes_local=isa_data['extra_includes_local'] + [os.path.join('include', os.path.basename(filename_decl))],
                                extra_includes_global=isa_data['extra_includes_global'],
                            )
                            filename_spec = os.path.join(
                                '../generated',
                                isa,
                                'lib',
                                '{}_{}_spec_{}_{}_{}_{}_{}.cpp'.format(
                                    kernel,
                                    isa,
                                    fptype,
                                    veclen,
                                    soalen,
                                    'compress18' if compress12 == 'true' else 'compress12',
                                    tbc_ts,
                                )
                            )
                            with open(filename_spec, 'w') as f:
                                f.write(rendered)


def _parse_args():
    '''
    Parses the command line arguments.

    :return: Namespace with arguments.
    :rtype: Namespace
    '''
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('isa', nargs='+')
    options = parser.parse_args()

    return options


if __name__ == "__main__":
    main()
