class Dog:
    def __init__(self, name):
        self.name = name


dog_1 = Dog("Arfie")

dog_2 = Dog("Spike")

dog_list = [dog_1, dog_2]

chosen_dog = [item for item in dog_list if item.name == "Spike"]

print(chosen_dog)
print(chosen_dog[0].name)