#!/usr/bin/env python3
"""
Policy YAML parser for state-change policies.

Reads policy.yaml and converts to structured dict for code generation.
"""

import yaml
from pathlib import Path
from typing import Dict, Any, List


class PolicyParser:
    """Parse and validate policy YAML files."""

    def __init__(self, yaml_path: Path):
        self.yaml_path = yaml_path
        self.policy = None

    def parse(self) -> Dict[str, Any]:
        """Parse YAML file and return structured policy dict."""
        with open(self.yaml_path, 'r') as f:
            self.policy = yaml.safe_load(f)

        # Basic validation
        self._validate_required_fields()

        # Resolve parameter references
        self._resolve_references()

        return self.policy

    def _validate_required_fields(self):
        """Check that required fields are present."""
        required = ['metadata', 'parameters', 'rules']

        for field in required:
            if field not in self.policy:
                raise ValueError(f"Missing required field: {field}")

        # Check metadata
        if 'name' not in self.policy['metadata']:
            raise ValueError("metadata.name is required")

        # Check rules
        if not self.policy['rules']:
            raise ValueError("At least one rule is required")

        for i, rule in enumerate(self.policy['rules']):
            if 'name' not in rule:
                raise ValueError(f"rules[{i}].name is required")
            if 'trigger' not in rule:
                raise ValueError(f"rules[{i}].trigger is required")
            if 'action' not in rule:
                raise ValueError(f"rules[{i}].action is required")

    def _resolve_references(self):
        """Resolve parameter references like 'parameter.cutoff'."""
        params = self.policy.get('parameters', {})

        def resolve_value(value):
            """Recursively resolve parameter references."""
            if isinstance(value, str) and value.startswith('parameter.'):
                param_name = value.split('.', 1)[1]
                if param_name not in params:
                    raise ValueError(f"Undefined parameter reference: {value}")
                return params[param_name]
            elif isinstance(value, dict):
                return {k: resolve_value(v) for k, v in value.items()}
            elif isinstance(value, list):
                return [resolve_value(v) for v in value]
            else:
                return value

        # Resolve in rules
        for rule in self.policy['rules']:
            rule['trigger'] = resolve_value(rule['trigger'])
            rule['action'] = resolve_value(rule['action'])
            if 'cooldown' in rule:
                rule['cooldown'] = resolve_value(rule['cooldown'])

    def get_fix_name(self) -> str:
        """Generate fix name from policy and problem."""
        problem = self.policy['metadata'].get('problem', 'generic')
        name = self.policy['metadata']['name']
        # Convert to CamelCase
        name_parts = name.replace('_', ' ').title().replace(' ', '')
        return f"FixStateChange{problem.capitalize()}{name_parts}"

    def get_cpp_filename(self) -> str:
        """Generate C++ filename."""
        problem = self.policy['metadata'].get('problem', 'generic')
        name = self.policy['metadata']['name']
        return f"fix_state_change_{problem}_{name}"


def main():
    """Test the parser."""
    import sys
    if len(sys.argv) < 2:
        print("Usage: python policy_parser.py <policy.yaml>")
        sys.exit(1)

    parser = PolicyParser(Path(sys.argv[1]))
    policy = parser.parse()

    print("Parsed policy:")
    print(f"  Name: {policy['metadata']['name']}")
    print(f"  Fix name: {parser.get_fix_name()}")
    print(f"  Filename: {parser.get_cpp_filename()}")
    print(f"  Rules: {len(policy['rules'])}")

    for i, rule in enumerate(policy['rules']):
        print(f"    {i+1}. {rule['name']}")
        print(f"       Trigger: {rule['trigger']['condition']}")
        print(f"       Action: {rule['action']['change_type']}")


if __name__ == '__main__':
    main()
