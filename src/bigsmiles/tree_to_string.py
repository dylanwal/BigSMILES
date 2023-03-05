"""

Code for turing a python class into a terminal printable tree.

"""
from typing import Protocol


class TreeConfig:
    tree_symbol_options = {
        'ascii': ('|', '|-- ', '+-- '),
        'ascii-ex': ('\u2502', '\u251c\u2500\u2500', '\u2514\u2500\u2500'),
        'ascii-exr': ('\u2502', '\u251c\u2500\u2500', '\u2570\u2500\u2500'),
        'ascii-em': ('\u2551', '\u2560\u2550\u2550', '\u255a\u2550\u2550'),
        'ascii-emv': ('\u2551', '\u255f\u2500\u2500', '\u2559\u2500\u2500'),
        'ascii-emh': ('\u2502', '\u255e\u2550\u2550', '\u2558\u2550\u2550'),
    }
    symbols = 'ascii-ex'  # some computers may need to show different symbols
    tab = "    "
    space = " "

    @classmethod
    def vertical_line(cls) -> str:
        return cls.tree_symbol_options[cls.symbols][0]

    @classmethod
    def tee(cls) -> str:
        return cls.tree_symbol_options[cls.symbols][1]

    @classmethod
    def corner(cls) -> str:
        return cls.tree_symbol_options[cls.symbols][2]


class TreeNode(Protocol):
    _tree_print_repr: bool


class TreeIntermediateNode(TreeNode):
    nodes: list


def to_repr(obj) -> str:
    return repr(obj)


def to_str(obj) -> str:
    return str(obj)


def tree_to_string(root_node: TreeIntermediateNode, show_object_labels: bool = True, print_repr: bool = False) -> str:
    """
    Creates a string representation of a tree for printing.

    Parameters
    ----------
    root_node: TreeIntermediateNode
        Root node
    show_object_labels: bool
        show object labels
    print_repr: bool
        print str or print repr of objects
        True: repr() is called
        False: str() is called

    Returns
    -------
    text: str
        String with the structure of a tree.

    """
    text = ""

    # add root node
    if show_object_labels:
        text = type(root_node).__name__ + ": "

    text += str(root_node)

    # add main part of tree
    text += tree_to_string_loop(root_node.nodes, [], show_object_labels, to_repr if print_repr else to_str)

    return text


def tree_to_string_loop(
        nodes: list[TreeNode, TreeIntermediateNode],
        spacers: list[bool],
        show_object_label: bool,
        func: callable
) -> str:
    """
    Creates the main body of the tree

    Parameters
    ----------
    nodes: list[TreeNode | TreeIntermediateNode]
        nodes
    spacers: list[bool]
        list records weather a continuation line (TreeConfig.vertical_line) is needed in the tree
        True: No line
        False: Add line
    show_object_label: bool
        show object/class type.__name__
    func: callable
        calls str() or repr()

    Returns
    -------
    text: str

    """
    text = ""

    for node in nodes:
        text += create_row(node, spacers, show_object_label, func, node == nodes[-1])

        # if interior node, do recursion
        if hasattr(node, "nodes"):
            text += tree_to_string_loop(
                node.nodes,
                spacers + [node == nodes[-1]],
                show_object_label,
                func
            )

    return text


def create_row(
        node: TreeNode | TreeIntermediateNode,
        spacers: list[bool],
        show_object_label: bool,
        func: callable,
        last_node: bool
) -> str:
    """
    Creates a row in the tree string representation.

    Parameters
    ----------
    node: TreeNode | TreeIntermediateNode
        node
    spacers: list[bool]
        list records weather a continuation line (TreeConfig.vertical_line) is needed in the tree
        True: No line
        False: Add line
    show_object_label: bool
        show object/class type.__name__
    func: callable
        calls str() or repr()
    last_node: bool
        last node in branch --> switch tee symbol to corner

    Returns
    -------
    text: str

    """
    line = "\n"

    # add spacer
    for i, spacer in enumerate(spacers):
        if not spacer:
            line += TreeConfig.vertical_line()
        line += TreeConfig.tab

    # add tree connection symbol
    if last_node:
        line += TreeConfig.corner() + TreeConfig.space
    else:
        line += TreeConfig.tee() + TreeConfig.space

    # add class label
    if show_object_label:
        line += type(node).__name__ + ": "

    # add node symbol
    if node._tree_print_repr:
        line += func(node)
    else:
        line += str(node)

    return line
