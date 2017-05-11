from vorbilt.bilayer_analyzer.bilayer_analyzer import print_valid_analyses

def test_print_valid_analyses():
    print("testing call to print_valid_analyses function...")
    print_valid_analyses()
    print(" ")
    print("without settings...")
    print_valid_analyses(show_settings=False)
    return

if __name__ == '__main__':
    test_print_valid_analyses()
