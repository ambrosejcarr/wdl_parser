import pyparsing as pp

_variable_chars = pp.alphanums + '._[]?'  # supported syntax in WDL
_value_chars = _variable_chars + '"'


def wdl_comment():
    return pp.Literal("#") + pp.restOfLine


def suppressed_keyword(s):
    return pp.Suppress(pp.Keyword(s))


def suppressed_literal(s):
    return pp.Suppress(pp.Literal(s))


def import_statement():
    return pp.Group(
        suppressed_keyword('import') +
        pp.QuotedString('"').setResultsName('imported_wdl') +
        pp.Optional(
            suppressed_keyword('as') +
            pp.Word(_variable_chars).setResultsName('imported_as_name')
        )
    )


def positional_argument():
    return pp.Group(
        pp.NotAny(pp.Keyword('call')) +
        pp.Word(_variable_chars).setResultsName('variable_type') +
        pp.Word(_variable_chars).setResultsName('variable_name')
    )


def keyword_argument():
    return pp.Group(
        pp.NotAny(pp.Keyword('call')) +
        pp.Regex('[a-zA-Z0-9_\[\]-]*\?').setResultsName('variable_type') +
        pp.Word(_variable_chars).setResultsName('variable_name') +
        pp.Optional(
            suppressed_literal('=') +
            pp.Word(_variable_chars).setResultsName('default_value')
        )
    )


def variable_definitions():
    return pp.ZeroOrMore(
        keyword_argument() |
        positional_argument()
    ).setResultsName('variable_definitions')


def call_input_assignment():
    return pp.Optional(
        suppressed_keyword('input:') +
        pp.OneOrMore(pp.Group(
            pp.Word(_variable_chars).setResultsName('variable_name') +
            suppressed_literal('=') +
            pp.Word(_variable_chars).setResultsName('value') +
            pp.Suppress(pp.Optional(pp.Literal(',')))))
    ).setResultsName('call_inputs')


def workflow_calls():
    return pp.Group(pp.OneOrMore(pp.Group(
        suppressed_keyword('call') +
        pp.Word(_variable_chars).setResultsName('task_name') +
        suppressed_literal('{') +
        call_input_assignment() +
        suppressed_literal('}')
    ))).setResultsName('calls')


def output_definition():
    return pp.OneOrMore(pp.Group(
        pp.Word(_variable_chars).setResultsName('variable_type') +
        pp.Word(_variable_chars).setResultsName('variable_name') +
        suppressed_literal('=') +
        pp.Word(_variable_chars + '"').setResultsName('variable_value')
    ))


def outputs():
    return (
        suppressed_keyword('output') +
        suppressed_literal('{') +
        output_definition() +
        suppressed_literal('}')
    ).setResultsName('outputs')


def workflow():
    return pp.Group(
        suppressed_keyword('workflow') +
        pp.Word(_variable_chars).setResultsName('workflow_name') +
        suppressed_literal('{') +
        variable_definitions() +
        workflow_calls() +
        outputs() +
        suppressed_literal('}')
    )


# grabs everything inside '{ }', aware of internal braces
def command():
    return (
        suppressed_literal('command') +
        pp.Group(pp.nestedExpr('{', '}')).setResultsName('command')
    )


def docker():
    return (
        suppressed_keyword('docker') +
        suppressed_literal(':') +
        pp.QuotedString('"').setResultsName('docker')
    )


def memory():
    return (
        suppressed_keyword('memory') +
        suppressed_literal(':') +
        pp.QuotedString('"').setResultsName('memory')
    )


def disks():
    return pp.Group(
        suppressed_keyword('disks') +
        pp.Suppress(':') +
        pp.Suppress('"') +
        pp.Word(pp.alphanums + '_/:-').setResultsName('disk_location') +
        pp.Word(pp.nums).setResultsName('size') +
        pp.Word(pp.alphas).setResultsName('disk_type') +
        pp.Suppress('"')
    ).setResultsName('disks')


def runtime():
    return (
        suppressed_literal('runtime') +
        suppressed_literal('{') +
        pp.Optional(docker()) +
        pp.Optional(memory()) +
        pp.Optional(disks()) +
        suppressed_literal('}')
    ).setResultsName('runtime')


def task():
    return pp.Group(
        suppressed_keyword('task') +
        pp.Word(_variable_chars).setResultsName('task_name') +
        suppressed_literal('{') +
        pp.Optional(variable_definitions()) +
        command() +
        pp.Optional(outputs()) +
        pp.Optional(runtime()) +
        suppressed_literal('}')
    )


def wdl():
    return pp.OneOrMore(
        task() |
        import_statement() |
        workflow()
    )




