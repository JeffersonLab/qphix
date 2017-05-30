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

    tf_list = [''.join(tfs)
               for tfs in itertools.product(*([('t', 's')] * 4))]
    true_false_list = [', '.join(tfs)
                       for tfs in itertools.product(*([('true', 'false')] * 4))]

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

        os.makedirs(os.path.join('..', 'generated', isa), exist_ok=True)

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
            print('Working on kernel `{}` for ISA `{}` …'.format(kernel, isa))
            defines = [
                (fptype, fptype_data['veclen'], fptype_data['soalens'])
                for fptype, fptype_data in sorted(isa_data['fptypes'].items())]

            # Generate the complete specialization.
            rendered = complete_specialization.render(
                generated_warning=generated_warning,
                ISA=isa,
                kernel_base=kernel,
                kernel_pattern=kernel_pattern,
                defines=defines,
                extra_includes_local=isa_data['extra_includes_local'],
                extra_includes_global=isa_data['extra_includes_global'],
                true_false_list=list(zip(true_false_list, tf_list)),
            )
            filename = os.path.join('../generated', isa, '{}_{}_complete_specialization.h'.format(kernel, isa))
            with open(filename, 'w') as f:
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
