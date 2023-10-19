# Creating New Elements in QuantumDotDesigner

QuantumDotDesigner is a versatile framework designed for the creation and manipulation of quantum dot devices. One of the key features of this framework is its modular design, allowing users to add new elements, each representing different components or functionalities within the quantum dot environment. This guide will walk you through the process of extending QuantumDotDesigner's capabilities by adding new elements.


## Understanding the Modular Design

The system is built around a core set of classes, each representing a different level of abstraction:

- `ElementBase`: This is the most basic type of object, containing properties and methods common to all elements.

- `Element`: This class extends `ElementBase`, adding further functionalities and properties.

- `Specific Elements (e.g., ClavierGate, Ohmic)`: These are subclasses of `Element`, each defining unique behaviour and properties for specific types of quantum dot components.

This hierarchy allows for easy expansion. To introduce a new type of element, you need to create a new class that inherits from `Element` and defines the unique characteristics and behaviour of the new component.

## Detailed Example: Adding a CustomGate Element

Suppose we want to create a `CustomGate` with a specific polygon shape. This element could represent a custom-designed gate in your quantum dot device. Below is a step-by-step guide on how to implement and integrate this new element.

### 1. Implementing the CustomGate Class

First, we define the `CustomGate` class, inheriting properties and methods from the `Element` class. We will override the `build` method to specify the geometry of our custom gate. Optionally, we can also generate the vertices from the newly created attributes in a separate method.

```python
from QuantumDotDesigner.base import Element
from QuantumDotDesigner.BaseCollection import BaseCollection
import gdstk

class CustomGate(Element):
    def __init__(self, name, collection: BaseCollection, vertices):
        super().__init__(name, collection)
        self.attr_1 = None
        self.attr_2 = None
        self.attr_3 = None
        # ... and more attributes if desired

	def get_vertices(self):
        """
        method to create polygon vertices from newly defined attributes 'self.attr_1',
        'self.attr_2', 'self.attr_3', ...
        """
        # code to calculate vertices
        # ...
        return vertices

    def build(self):
        if self.layer is None or not isinstance(self.layer, Layer):
            raise ValueError("A valid Layer must be assigned before building the element.")
		vertices = self.get.vertices()
        polygon = gdstk.Polygon(vertices, layer=self.layer)
        cell = gdstk.Cell(self.name)
        cell.add(polygon)
        self.elements[self.name]['vertices'] = polygon.points
        self.elements[self.name]['positions'] = [[0, 0]]
        self.elements[self.name]['layer'] = self.layer
        self.elements[self.name]['layer_stage'] = self._layer_stage
        self.cell = cell
        self._set_built(True)
```

### 2. Using the CustomGate Class

After defining our `CustomGate`, we need to create an instance of it and integrate it into our quantum dot design.

```python

from QuantumDotDesigner.BaseCollection import BaseCollection

# other imports...

collection = BaseCollection()


custom_gate = CustomGate(name="my_custom_gate", collection=collection)

# assign the attributes to your created instance
custom_gate.attr_1 = 10
custom_gate.attr_2 = True
custom_gate.attr_3 = 'type_x'

custom_gate.build()
```

### Conclusion

By following these steps, you have created a new type of element with a specific geometric shape and integrated it into your quantum dot design. This process showcases the flexibility and extensibility of the QuantumDotDesigner framework, allowing for the incorporation of custom components and functionalities.

---
