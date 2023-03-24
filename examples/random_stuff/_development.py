
import bigsmiles.constructors.tokenizer


def generate_list_of_tokens():
    text = "C%123"
    results = bigsmiles.constructors.tokenizer.tokenize(text)

    output = "[" + ", ".join([str(token.kind) for token in results]) + "]"
    print(output)


if __name__ == "__main__":
    generate_list_of_tokens()
