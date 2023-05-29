from super_hash import *

class THingy():
    pass


exec("""
def thing():
    pass
""")

def static_nice_function():
    pass
    
print(function_hashers.smart(static_nice_function,debug=True))
print(function_hashers.smart(thing,debug=True))
print(function_hashers.smart(THingy,debug=True))
print(function_hashers.shallow(static_nice_function))
print(function_hashers.shallow(thing))

# from super_hash import super_hash as hash
value = {
    frozenset({
        frozenset({
            "key-deep-deep": 10
        }.items()): "key-deep",
    }.items()): "first_value",
    "second_value": [
        {"a": 10},
    ]
}
print(super_hash(value))