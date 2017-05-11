from vorbilt.bilayer_analyzer.bilayer_analyzer import valid_analyses

def test_valid_analyses():
    print("testing call to valid_analyses function...")
    valid_analyses()
    print(" ")
    print("without settings...")
    valid_analyses(show_settings=False)
    return

if __name__ == '__main__':
    test_valid_analyses()
