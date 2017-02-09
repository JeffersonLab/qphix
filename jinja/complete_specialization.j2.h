{# Copyright © 2017 Martin Ueding <dev@martin-ueding.de> #}

{% import 'jinja/common.j2.h' as common %}

#include "qphix/geometry.h"
#include "qphix/avx/avx_utils.h"

{{ common.head() }}
{{ common.general(kernel) }}

{% for fptype, veclen, soalens in defines %}
{{ common.specialization(kernel, isa, fptype, veclen, soalens) }}
{% endfor %}
