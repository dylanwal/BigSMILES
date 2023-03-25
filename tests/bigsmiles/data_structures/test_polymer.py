import pytest



# def test_extract_text_between_pipe_no_pipes():
#     string = "This string has no pipes"
#     removed_text, remaining_text = extract_text_between_pipe(string)
#     assert removed_text == []
#     assert remaining_text == [string]
#
#
# def test_extract_text_between_pipe_single_pipe():
#     string = "This string has |one pipe|"
#     removed_text, remaining_text = extract_text_between_pipe(string)
#     assert removed_text == ["one pipe"]
#     assert remaining_text == ["This string has ", ""]
#
#
# def test_extract_text_between_pipe_multiple_pipes():
#     string = "|First| string with |multiple| pipes |in between|"
#     removed_text, remaining_text = extract_text_between_pipe(string)
#     assert removed_text == ["First", "multiple", "in between"]
#     assert remaining_text == ["", " string with ", " pipes ", ""]
#
#
# def test_extract_text_between_pipe_consecutive_pipes():
#     string = "This string || has consecutive pipes ||"
#     removed_text, remaining_text = extract_text_between_pipe(string)
#     assert removed_text == ["", " has consecutive pipes ", ""]
#     assert remaining_text == ["This string ", " ", ""]
#
#
# def test_extract_text_between_pipe_only_pipes():
#     string = "||"
#     removed_text, remaining_text = extract_text_between_pipe(string)
#     assert removed_text == ["", ""]
#     assert remaining_text == ["", ""]
