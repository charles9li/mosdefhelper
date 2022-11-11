__all__ = ['PolyAlkylAcrylateUA']
from ast import literal_eval
import string

import mbuild as mb

import mbuildhelper.lib.molecules as mbh_lib_molecules


_OPERATOR_PRECEDENCE = {'+': 2,
                        '*': 3}
_OPERATORS = list(_OPERATOR_PRECEDENCE.keys())
_DELIMITERS = ['(', ')']


def _tokenize_expr(expr):
    token_list = expr.split()
    for char in _OPERATORS + _DELIMITERS:
        temp_list = []
        for token in token_list:
            token = [t.strip() for t in token.split(char)]
            token_split = [char] * (len(token) * 2 - 1)
            token_split[0::2] = token
            temp_list += token_split
        token_list = temp_list
    return [t for t in token_list if t]


def _create_stack(token_list):
    output = []
    stack = []
    for token in token_list:
        if token in _OPERATORS:
            while len(stack) > 0 and stack[-1] != '(' and _OPERATOR_PRECEDENCE[stack[-1]] >= _OPERATOR_PRECEDENCE[token]:
                output.append(stack.pop())
            stack.append(token)
        elif token == '(':
            stack.append(token)
        elif token == ')':
            while stack[-1] != '(':
                output.append(stack.pop())
            stack.pop()
        else:
            output.append(token)
    while len(stack) > 0:
        output.append(stack.pop())
    return output


def _evaluate_output(output):
    stack = []
    for item in output:
        if item in _OPERATORS:
            op2 = stack.pop()
            op1 = stack.pop()
            if item == '+':
                stack.append(op1 + op2)
            elif item == '*':
                stack.append(op1 * op2)
        else:
            try:
                stack.append(int(item))
            except ValueError:
                stack.append([item])
    if len(stack) != 1:
        raise ValueError("Malformed expression.")
    return stack.pop()


def _expression_to_sequence(expr):
    return _evaluate_output(_create_stack(_tokenize_expr(expr)))


def _token_to_monomer(monomer_token, cap_front=False, cap_end=False):
    if monomer_token.startswith("A"):
        alkyl_tail_length = literal_eval(monomer_token[1:])
        methyl = False
    elif monomer_token.startswith("mA"):
        alkyl_tail_length = literal_eval(monomer_token[2:])
        methyl = True
    else:
        raise ValueError(f"invalid monomer token: {monomer_token}")
    return mbh_lib_molecules.AlkylAcrylateUA(n=alkyl_tail_length, methyl=methyl, cap_front=cap_front, cap_end=cap_end)


_CHARACTERS = string.ascii_letters


class PolyAlkylAcrylateUA(mb.Compound):

    def __init__(self, expression, **kwargs):
        super(PolyAlkylAcrylateUA, self).__init__(**kwargs)

        # parse expression to list of monomers
        monomer_sequence = _expression_to_sequence(expression)

        # single molecule case
        if len(monomer_sequence) == 1:
            self.add(_token_to_monomer(monomer_sequence[0], cap_front=True, cap_end=True), label="monomer[$]")
        # polymer case
        else:
            # get the front and end monomers from sequence
            front_monomer = _token_to_monomer(monomer_sequence.pop(0), cap_front=True)
            end_monomer = _token_to_monomer(monomer_sequence.pop(-1), cap_end=True)

            # add monomers
            last_part = front_monomer
            self.add(last_part, label="monomer[$]")
            for m in monomer_sequence:
                this_part = _token_to_monomer(m)
                self.add(this_part, label="monomer[$]")
                mb.force_overlap(
                    move_this=this_part,
                    from_positions=this_part.labels['up'],
                    to_positions=last_part.labels['down']
                )
                last_part = this_part
            self.add(end_monomer, label="monomer[$]")
            mb.force_overlap(
                move_this=end_monomer,
                from_positions=end_monomer.labels['up'],
                to_positions=last_part.labels['down']
            )


if __name__ == '__main__':
    pass
