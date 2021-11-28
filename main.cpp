#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <basix/finite-element.h>

#include <memory>

using namespace dealii;

template <int dim, int spacedim = dim>
class BasixFE : public FiniteElement<dim, spacedim>
{
private:
  static FiniteElementData<dim>
  generate_finite_element_data(const std::unique_ptr<basix::FiniteElement> &fe)
  {
    std::vector<unsigned int> dpo(dim + 1);

    for (unsigned int i = 0; i <= dim; ++i)
      dpo[i] = fe->num_entity_dofs()[i][0 /*TODO*/];

    ReferenceCell reference_cell;

    switch (fe->cell_type())
      {
        case basix::cell::type::triangle:
          reference_cell = ReferenceCells::Triangle;
          break;
        default:
          Assert(false, ExcNotImplemented());
      }

    const unsigned int n_components = 1; // TODO

    const unsigned int degree = fe->degree();

    typename FiniteElementData<dim>::Conformity conformity;

    if (fe->discontinuous())
      conformity = FiniteElementData<dim>::Conformity::L2;
    else
      conformity = FiniteElementData<dim>::Conformity::H1;

    return FiniteElementData<dim>(
      dpo, reference_cell, n_components, degree, conformity);
  }


public:
  BasixFE(const std::unique_ptr<basix::FiniteElement> &fe)
    : FiniteElement<dim, spacedim>(
        generate_finite_element_data(fe),
        std::vector<bool>(fe->dim()),
        std::vector<ComponentMask>(fe->dim(), ComponentMask(1 /*TODO*/, true)))
    , fe(fe)
  {
    if (fe->dof_transformations_are_identity())
      {
        // nothing to do
      }
    else if (fe->dof_transformations_are_permutations())
      {
        AssertDimension(dim, 2);

        const auto n_vertices = this->reference_cell().n_vertices();
        const auto n_lines    = this->reference_cell().n_lines();

        orienation_table.resize(n_lines);

        unsigned int counter = n_vertices * this->n_dofs_per_vertex();

        const auto &trafo = fe->base_transformations();

        for (unsigned int l = 0; l < n_lines; ++l)
          {
            // TODO: deal.II assumes that all faces have the same number of
            // DoFs
            const unsigned int n_dofs = this->n_dofs_per_line();

            orienation_table[l].resize(n_dofs);

            for (unsigned int i = 0; i < n_dofs; ++i)
              for (unsigned int j = 0; j < n_dofs; ++j)
                {
                  if (trafo(l, counter + j, counter + i))
                    {
                      orienation_table[l][i] = j;
                      break;
                    }
                }

            counter += n_dofs;
          }
      }
    else
      {
        Assert(false, ExcNotImplemented());
      }
  }

  std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override
  {
    return std::make_unique<BasixFE<dim, spacedim>>(*this);
  }

  std::string
  get_name() const override
  {
    std::ostringstream namebuf;
    namebuf << "BasixFE<" << dim << ">(" << this->degree << ")";

    return namebuf.str();
  }

  virtual unsigned int
  adjust_quad_dof_index_for_face_orientation(
    const unsigned int index,
    const unsigned int face_no,
    const bool         face_orientation,
    const bool         face_flip,
    const bool         face_rotation) const override
  {
    if (fe->dof_transformations_are_identity())
      return index;

    Assert(false, ExcNotImplemented());

    (void)index;
    (void)face_no;
    (void)face_orientation;
    (void)face_flip;
    (void)face_rotation;
  }

  virtual unsigned int
  adjust_line_dof_index_for_line_orientation(
    const unsigned int index,
    const bool         line_orientation) const override
  {
    if (fe->dof_transformations_are_identity())
      return index;

    if (line_orientation == true)
      return index;

    // TODO: deal.II assumes that all faces have the same number of DoFs
    const unsigned int line_no = 0;

    return orienation_table[line_no][index];
  }

  UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override
  {
    Assert(false, ExcNotImplemented());

    (void)update_flags;

    return {};
  }

  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim> &       quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override
  {
    Assert(false, ExcNotImplemented());

    (void)update_flags;
    (void)mapping;
    (void)quadrature;
    (void)output_data;

    return {};
  }

  void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const
  {
    Assert(false, ExcNotImplemented());

    (void)cell;
    (void)cell_similarity;
    (void)quadrature;
    (void)mapping;
    (void)mapping_internal;
    (void)mapping_data;
    (void)fe_internal;
    (void)output_data;
  }

  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override
  {
    Assert(false, ExcNotImplemented());

    (void)cell;
    (void)face_no;
    (void)sub_no;
    (void)quadrature;
    (void)mapping;
    (void)mapping_internal;
    (void)mapping_data;
    (void)fe_internal;
    (void)output_data;
  }

private:
  const std::unique_ptr<basix::FiniteElement> &fe;
  std::vector<std::vector<unsigned int>>       orienation_table;
};

/**
 * /home/munch/sw-basix/cmake-3.22.0-install/bin/cmake ..
 * -DBasix_DIR=/home/munch/sw-basix/basix-install/lib64/cmake/basix/
 * -Dxtl_DIR=/home/munch/sw-basix/basix/_deps/xtl-build
 * -Dxtensor_DIR=/home/munch/sw-basix/basix/_deps/xtensor-build
 * -DDEAL_II_DIR=../dealii-build
 */
int
main()
{
  auto _element = std::make_unique<basix::FiniteElement>(
    basix::create_element(basix::element::family::P,
                          basix::cell::type::triangle,
                          3,
                          basix::element::lagrange_variant::equispaced,
                          false));

  std::cout << _element->dim() << std::endl;

  for (int j = 0; j < _element->dim(); ++j)
    {
      for (unsigned int i = 0; i < 2; ++i)
        std::cout << _element->points()(j, i) << " ";
      std::cout << std::endl;
    }

  const auto &trafo = _element->base_transformations();

  std::cout << trafo.shape()[0] << std::endl;
  std::cout << trafo.shape()[1] << std::endl;
  std::cout << trafo.shape()[2] << std::endl;

  for (unsigned int i = 0; i < trafo.shape()[0]; ++i)
    {
      for (unsigned int j = 0; j < trafo.shape()[1]; ++j)
        {
          for (unsigned int k = 0; k < trafo.shape()[2]; ++k)
            std::cout << trafo(i, j, k) << " ";
          std::cout << std::endl;
        }
      std::cout << std::endl;
    }

  {
    const unsigned int dim = 2;

    Triangulation<dim> tria;

    GridGenerator::subdivided_hyper_rectangle_with_simplices(tria,
                                                             {1, 1},
                                                             {0.0, 0.0},
                                                             {1.0, 1.0});

    DoFHandler<dim> dof_handler(tria);

    BasixFE<dim> fe(_element);

    dof_handler.distribute_dofs(fe);

    std::vector<types::global_dof_index> dof_indices;

    for (const auto cell : dof_handler.active_cell_iterators())
      {
        dof_indices.resize(cell->get_fe().n_dofs_per_cell());

        cell->get_dof_indices(dof_indices);

        for (const auto i : dof_indices)
          std::cout << std::setw(3) << i << " ";
        std::cout << std::endl;
      }
  }

  return 0;
}
