import pyparsing as pp

# todo maybe need comments throughout

_valid_chars = pp.alphanums + '._[]?'  # supported syntax in WDL

# start by parsing workflows
# todo need to parse default arguments here (look for =, if found, parse next)
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

_command = (
    pp.Suppress(pp.Literal('command')) +
    pp.nestedExpr('{', '}')  # grabs everything inside '{ }', aware of internal braces
).setResultsName('command')

_docker = (
    pp.Suppress(pp.Keyword('docker')) +
    pp.Suppress(pp.Literal(':')) +
    pp.Word(pp.alphanums + '_/:')
).setResultsName('docker')

_memory = (
    pp.Suppress(pp.Keyword('memory')) +
    pp.Suppress(pp.Literal(':')) +
    pp.Word(pp.alphanums + '_/:')
).setResultsName('memory')

_disk = pp.Group(
    pp.Suppress(pp.Keyword('disks')) +
    pp.Suppress(pp.Literal(':')) +
    pp.Word(pp.alphanums + '_/:').setResultsName('disk_location') +
    pp.Word(pp.nums).setResultsName('size') +
    pp.Word(pp.alphas).setResultsName('disk_type')
).setResultsName('disks')

_runtime = (
    pp.Suppress(pp.Literal('runtime')) +
    pp.Suppress(pp.Literal('{')) +
    pp.Optional(_docker) +
    pp.Optional(_memory) +
    pp.Optional(_disk) +
    pp.Suppress(pp.Literal('}'))
).setResultsName('runtime')

_task = (
    pp.Suppress(pp.Keyword('task')) +
    pp.Word(_valid_chars).setResultsName('task_name') +
    pp.Suppress(pp.Literal('{')) +
    pp.Optional(_variable_definitions) +
    _command +
    pp.Optional(_runtime) +
    pp.Optional(_outputs) +
    pp.Suppress(pp.Literal('}'))
)

_imports = NotImplemented































