import os
from . import parse

# todo might be able to switch to getattr/setattr for code maintainability


class ImportedWDL:

    def __init__(self, parsed_import):
        self._name = parsed_import['imported_wdl']
        if 'imported_as_name' in parsed_import:
            self._alias = parsed_import['imported_as_name']


class Input:

    def __init__(self, parsed_input):
        if 'default_value' in parsed_input:
            self._default_value = parsed_input['default_value']
        else:
            self._default_value = None
        self._type = parsed_input['variable_type']
        self._name = parsed_input['variable_name']


class Output:

    def __init__(self, parsed_output):
        self._type = parsed_output['variable_type']
        self._name = parsed_output['variable_name']
        self._value = parsed_output['variable_value']


class Call:

    def __init__(self, parsed_call):
        self._name = parsed_call['task_name']
        self._inputs = [Input(v) for v in parsed_call['call_inputs']]


class Workflow:

    def __init__(self, parsed_workflow):
        self._name = parsed_workflow['workflow_name']
        self._inputs = [Input(v) for v in parsed_workflow['variable_definitions']]
        self._calls = [Call(v) for v in parsed_workflow['calls']]
        self._outputs = [Output(v) for v in parsed_workflow['output']]


class Command:

    def __init__(self, parsed_command):
        self._command = parsed_command


class Disks:

    def __init__(self, parsed_disks):
        self._location = parsed_disks['disk_location']
        self._size = parsed_disks['size']
        self._type = parsed_disks['disk_type']


class Memory:

    def __init__(self, memory):
        self._value = memory


class Docker:

    def __init__(self, docker_name):
        self._name = docker_name


class Runtime:

    def __init__(self, parsed_runtime):
        try:
            self._docker = Docker(parsed_runtime['docker'])
        except KeyError:
            self._docker = None
        try:
            self._memory = Memory(parsed_runtime['memory'])
        except KeyError:
            self._memory = None
        try:
            self._disks = Disks(parsed_runtime['disks'])
        except KeyError:
            self._disks = None


class Task:

    def __init__(self, parsed_task):
        self._name = parsed_task['task_name']
        self._inputs = [Input(v) for v in parsed_task['variable_definitions']]
        self._command = parsed_task['command']

        if 'outputs' in parsed_task:
            self._outputs = [Output(v) for v in parsed_task['outputs']]
        else:
            self._outputs = None

        if 'runtime' in parsed_task:
            self._runtime = Runtime(parsed_task['runtime'])
        else:
            self._runtime = None


class WDL:

    def __init__(self, wdl):

        if os.path.isfile(wdl):
            with open(wdl, 'r') as f:
                data = f.read()
        else:
            data = wdl

        self._tasks = []
        self._workflows = []
        self._imports = []

        parsed = parse.wdl().ignore(parse.wdl_comment()).parseString(data)

        for block in parsed:
            if 'task_name' in block:
                self._tasks.append(Task(block))
            elif 'workflow_name' in block:
                self._workflows.append(Workflow(block))
            elif 'imported_wdl' in block:
                self._imports.append(ImportedWDL(block))

    @property
    def tasks(self):
        return self._tasks

    @property
    def workflows(self):
        return self._workflows

    @property
    def imports(self):
        return self._imports
