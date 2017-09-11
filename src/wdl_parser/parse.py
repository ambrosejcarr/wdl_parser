import pyparsing as pp

_valid_chars = pp.alphanums + '._[]'  # supported syntax in WDL

# start by parsing workflows
_variable_definitions = pp.ZeroOrMore(pp.Group(
    pp.NotAny(pp.Keyword('call')) +
    pp.Word(_valid_chars).setResultsName('type') +
    pp.Word(_valid_chars).setResultsName('variable_name'))
).setResultsName('variable_definitions')

_call_assigned_input = (
    pp.Suppress(pp.Keyword('input:')) +
    pp.ZeroOrMore(pp.Group(
        pp.Word(_valid_chars).setResultsName('variable_name') +
        pp.Suppress(pp.Literal('=')) +
        pp.Word(_valid_chars).setResultsName('value') +
        pp.Suppress(pp.Optional(pp.Literal(',')))))
).setResultsName('inputs')

_calls = pp.Group(pp.OneOrMore(pp.Group(
    pp.Suppress(pp.Keyword('call')) +
    pp.Word(_valid_chars).setResultsName('task_name') +
    pp.Suppress(pp.Literal('{')) +
    _call_assigned_input +
    pp.Suppress(pp.Literal('}'))
))).setResultsName('calls')

_output_object_definition = pp.OneOrMore(pp.Group(
    pp.Word(_valid_chars).setResultsName('variable_type') +
    pp.Word(_valid_chars).setResultsName('variable_name') +
    pp.Suppress(pp.Literal('=')) +
    pp.Word(_valid_chars).setResultsName('variable_value')
))

_outputs = (
    pp.Suppress(pp.Keyword('output')) +
    pp.Suppress(pp.Literal('{')) +
    _output_object_definition +
    pp.Suppress(pp.Literal('}'))
).setResultsName('outputs')

_workflow = (
    pp.Suppress(pp.Keyword('workflow')) +
    pp.Word(_valid_chars).setResultsName('workflow_name') +
    pp.Suppress(pp.Literal('{')) +
    _variable_definitions +
    _calls +
    _outputs +
    pp.Suppress(pp.Literal('}'))
)
